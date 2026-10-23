"""Loopback-only local research application; serves only the public directory."""
from http.server import ThreadingHTTPServer, SimpleHTTPRequestHandler
from pathlib import Path
from urllib.parse import urlparse
import argparse,json,traceback
from physics import waveform,fit_event,DEFAULT

ROOT=Path(__file__).resolve().parent

class Handler(SimpleHTTPRequestHandler):
    def __init__(self,*args,**kwargs):super().__init__(*args,directory=str(ROOT/'public'),**kwargs)
    def send_json(self,obj,status=200):
        body=json.dumps(obj,allow_nan=False,separators=(',',':')).encode()
        self.send_response(status);self.send_header('Content-Type','application/json; charset=utf-8');self.send_header('Content-Length',str(len(body)));self.send_header('Cache-Control','no-store');self.end_headers();self.wfile.write(body)
    def do_GET(self):
        path=urlparse(self.path).path
        if path=='/api/health':return self.send_json({'ok':True})
        if path=='/api/event':
            f=ROOT/'data/event.json'
            return self.send_json(json.loads(f.read_text()) if f.exists() else {'error':'Run prepare_event.py first'},200 if f.exists() else 503)
        if path.startswith('/api/'):return self.send_json({'error':'Unknown route'},404)
        return super().do_GET()
    def do_POST(self):
        try:
            origin=self.headers.get('Origin')
            if origin and origin!=f'http://{self.headers.get("Host")}':return self.send_json({'error':'Origin rejected'},403)
            size=int(self.headers.get('Content-Length','0'))
            if not 0<size<8192:raise ValueError('Invalid request size')
            obj=json.loads(self.rfile.read(size))
            if self.path=='/api/waveform':return self.send_json(waveform(obj))
            if self.path=='/api/fit':return self.send_json(fit_event(obj.get('parameters',{}),obj.get('detector','H1')))
            return self.send_json({'error':'Unknown route'},404)
        except (ValueError,KeyError,TypeError) as exc:self.send_json({'error':str(exc)},400)
        except Exception:
            traceback.print_exc();self.send_json({'error':'Calculation failed; see local server log'},500)
    def log_message(self,fmt,*args):
        if '/api/health' not in str(args):print(fmt%args,flush=True)

if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--port',type=int,default=8765);args=parser.parse_args()
    print('Warming waveform model…',flush=True)
    waveform(DEFAULT)
    print(f'Local: http://127.0.0.1:{args.port}',flush=True)
    ThreadingHTTPServer(('127.0.0.1',args.port),Handler).serve_forever()
