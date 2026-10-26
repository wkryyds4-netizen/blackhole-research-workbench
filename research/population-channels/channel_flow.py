"""Original adapter for the published CC-BY4 channel emulator tensors.

The package supplies the network mathematics; only numeric state tensors and
published coordinate mappings are read. No AMA​ZE source module is executed.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
import torch
from scipy.special import expit,logit
from glasflow import CouplingNSF

BASE=Path(__file__).resolve().parent
PRODUCTS=BASE/'sources/inspected-products'
CHANNELS=['CE','CHE','GC','NSC','SMT']


def sha(path):
    with path.open('rb') as stream:return hashlib.file_digest(stream,'sha256').hexdigest()


class ChannelFlow:
    def __init__(self,channel):
        assert channel in CHANNELS
        report=json.loads((BASE/'product-inspection.json').read_text())
        for name in ['flowconfig.json','flow_mappings.json']:
            assert json.loads((PRODUCTS/name).read_text())==report['json_metadata'][name]
        config=report['json_metadata']['flowconfig.json'][channel]
        mapping=report['json_metadata']['flow_mappings.json'][channel]
        self.channel=channel;self.context_size=2 if channel=='CE' else 1
        self.upper=np.array([mapping['Mup'],mapping['qup'],1.,mapping['zup']])
        self.scale=np.array([mapping['Mmax'],mapping['qmax'],1.,mapping['zmax']])
        path=PRODUCTS/(channel+'.pt')
        assert sha(path)==report['torch_archives'][channel+'.pt']['sha256']
        self.network=CouplingNSF(n_inputs=4,n_conditional_inputs=self.context_size,
            n_transforms=config['transforms'],n_neurons=config['neurons'],n_blocks_per_transform=2,
            batch_norm_between_transforms=True,linear_transform=None,num_bins=config['bins'])
        state=torch.load(path,map_location='cpu',weights_only=True)
        assert isinstance(state,dict) and all(isinstance(k,str) and isinstance(v,torch.Tensor) for k,v in state.items())
        assert all(torch.isfinite(v).all() for v in state.values() if torch.is_floating_point(v))
        loaded=self.network.load_state_dict(state,strict=True)
        assert not loaded.missing_keys and not loaded.unexpected_keys
        self.state_keys=len(state)
        self.network.double().eval()

    def to_network(self,theta):
        theta=np.asarray(theta,float);assert theta.ndim==2 and theta.shape[1]==4
        assert np.isfinite(theta).all()
        for i in [0,1,3]:assert np.all((theta[:,i]>0)&(theta[:,i]<self.upper[i]))
        assert np.all(abs(theta[:,2])<1)
        mapped=np.empty_like(theta)
        for i in [0,1,3]:mapped[:,i]=logit(theta[:,i]/self.upper[i])/self.scale[i]
        mapped[:,2]=np.arctanh(theta[:,2])
        return mapped

    def to_physical(self,mapped):
        mapped=np.asarray(mapped,float);theta=np.empty_like(mapped)
        for i in [0,1,3]:theta[:,i]=self.upper[i]*expit(mapped[:,i]*self.scale[i])
        theta[:,2]=np.tanh(mapped[:,2])
        return theta

    def log_coordinate_jacobian(self,theta):
        theta=np.asarray(theta,float)
        jac=np.empty_like(theta)
        for i in [0,1,3]:jac[:,i]=self.upper[i]/(theta[:,i]*(self.upper[i]-theta[:,i])*self.scale[i])
        jac[:,2]=1/(1-theta[:,2]**2)
        return np.log(abs(jac)).sum(axis=1)

    def context(self,raw_context,n):
        raw=np.asarray(raw_context,float)
        assert raw.shape==(self.context_size,) and np.isfinite(raw).all()
        return torch.as_tensor(np.broadcast_to(raw,(n,self.context_size)).copy(),dtype=torch.float64)

    def log_density(self,theta,raw_context):
        mapped=self.to_network(theta)
        with torch.no_grad():values=self.network.log_prob(torch.as_tensor(mapped),conditional=self.context(raw_context,len(mapped))).numpy()
        return values+self.log_coordinate_jacobian(theta)

    def from_latent(self,latent,raw_context):
        latent=np.asarray(latent,float)
        with torch.no_grad():mapped,_=self.network.inverse(torch.as_tensor(latent),conditional=self.context(raw_context,len(latent)))
        return self.to_physical(mapped.numpy())
