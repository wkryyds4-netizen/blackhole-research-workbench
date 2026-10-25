"""Strict bounded HTTP Range reader for read-only HDF metadata inspection."""
import io
import hashlib
import time
import json
from pathlib import Path
import requests

class RateLimitError(IOError):
    pass


class RangeReader(io.RawIOBase):
    def __init__(self,url,size,max_bytes=24*1024*1024,max_requests=256,block_size=65536,cache_key=None,min_interval=0):
        super().__init__();self.url=url;self.size=size;self.position=0
        self.max_bytes=max_bytes;self.max_requests=max_requests;self.block_size=block_size
        self.session=requests.Session();self.cache={};self.requests=0;self.network_bytes=0;self.ranges=[];self.identity=None
        self.cached_bytes=0
        self.min_interval=min_interval;self.last_request=0
        self.disk_cache=None if cache_key is None else Path(__file__).resolve().parent/'.http-cache'/hashlib.sha256((url+'|'+str(size)+'|'+cache_key).encode()).hexdigest()
        if self.disk_cache is not None:self.disk_cache.mkdir(parents=True,exist_ok=True)

    def readable(self):return True
    def seekable(self):return True
    def tell(self):return self.position

    def seek(self,offset,whence=io.SEEK_SET):
        position=offset if whence==io.SEEK_SET else self.position+offset if whence==io.SEEK_CUR else self.size+offset
        if position<0:raise ValueError('Negative seek')
        self.position=position;return position

    def block(self,start):
        if start in self.cache:return self.cache[start]
        end=min(self.size-1,start+self.block_size-1);length=end-start+1
        if sum(len(v) for v in self.cache.values())+length>self.max_bytes:raise IOError('Logical range byte budget exhausted')
        if self.disk_cache is not None:
            cache=self.disk_cache/f'{start}-{end}.bin';meta=cache.with_suffix('.json')
            if cache.exists() and meta.exists():
                content=cache.read_bytes();record=json.loads(meta.read_text())
                if len(content)!=length or hashlib.sha256(content).hexdigest()!=record['sha256']:raise IOError('Cached range failed checksum')
                self.cache[start]=content;self.ranges.append(record);self.cached_bytes+=length
                return content
        if self.requests>=self.max_requests or self.network_bytes+length>self.max_bytes:raise IOError('Metadata HTTP range budget exhausted')
        for attempt in range(3):
            self.requests+=1
            if self.requests>self.max_requests:raise IOError('Metadata request budget exhausted')
            try:
                delay=self.min_interval-(time.monotonic()-self.last_request)
                if delay>0:time.sleep(delay)
                self.last_request=time.monotonic()
                with self.session.get(self.url,headers={'Range':f'bytes={start}-{end}','Accept-Encoding':'identity'},stream=True,timeout=(10,25)) as response:
                    if response.status_code==429:raise RateLimitError('HTTP 429; Retry-After='+str(response.headers.get('Retry-After'))+'; stop new requests')
                    response.raise_for_status()
                    if response.status_code!=206 or response.headers.get('Content-Range')!=f'bytes {start}-{end}/{self.size}':raise IOError('Server did not honor exact range')
                    identity=response.headers.get('ETag') or response.headers.get('Last-Modified')
                    if self.identity is None:self.identity=identity
                    if identity!=self.identity:raise IOError('Remote object identity changed')
                    chunks=[];received=0
                    for data in response.iter_content(65536):
                        received+=len(data);self.network_bytes+=len(data)
                        if received>length or self.network_bytes>self.max_bytes:raise IOError('Range byte budget exceeded')
                        chunks.append(data)
                    content=b''.join(chunks)
                    if len(content)!=length:raise IOError('Truncated HTTP range')
                self.cache[start]=content
                record=dict(start=start,end=end,sha256=hashlib.sha256(content).hexdigest())
                self.ranges.append(record)
                if self.disk_cache is not None:
                    cache.write_bytes(content);meta.write_text(json.dumps(record),encoding='utf-8')
                return content
            except requests.RequestException:
                if attempt==2:raise
                time.sleep(.3*(attempt+1))

    def readinto(self,b):
        count=min(len(b),max(0,self.size-self.position));written=0
        while written<count:
            start=(self.position//self.block_size)*self.block_size;data=self.block(start)
            offset=self.position-start;n=min(count-written,len(data)-offset)
            b[written:written+n]=data[offset:offset+n];written+=n;self.position+=n
        return written

    def close(self):
        self.session.close();super().close()
