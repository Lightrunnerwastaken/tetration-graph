"""Independent research workspace; the upstream projects and caches are read only."""
import hashlib,json,sys,time,queue
from pathlib import Path
sys.dont_write_bytecode=True
WS=Path(__file__).absolute().parents[2]
HERE=Path(__file__).absolute().parent
PROJECT=Path(r"D:\Users\janis\Documents\Matura\tetration-graph")
GP_SRC=Path(r"D:\Users\janis\Documents\Tetration\gp-tetration\src")
sys.path.insert(0,str(HERE/'_backend') if (HERE/'_backend/fatou_backend').exists() else str(GP_SRC))
import mpmath as mp
from fatou_backend import FatouGP,find_default_gp_exe
from fatou_backend.worker import FatouGPWorker,WorkerDied

class LoggedWorker(FatouGPWorker):
    """Stream construction diagnostics; kill only this worker on its deadline."""
    def _read_until(self,sentinel,timeout):
        deadline=time.perf_counter()+timeout;collected=[];errors=[]
        while True:
            remaining=deadline-time.perf_counter()
            if remaining<=0:
                self.close()
                raise WorkerDied(f"timeout after {timeout}s; last output:\n"+"\n".join(collected[-30:]))
            try:line=self._queue.get(timeout=min(remaining,1))
            except queue.Empty:continue
            if line is None:
                self.close();raise WorkerDied("worker ended; last output:\n"+"\n".join(collected[-30:]))
            if line.strip()==sentinel:
                if errors:
                    self.close();raise RuntimeError("PARI/GP error:\n"+"\n".join(errors))
                return collected
            collected.append(line)
            if "o8kn:" in line or "hinit:" in line:print(line,flush=True)
            if "***" in line and "Warning" not in line:errors.append(line)

class LoggedGP(FatouGP):
    def _spawn_full(self,base):
        return LoggedWorker(self.gp_exe,self._init_lines(base),self.init_timeout,self.eval_timeout)
    def _run_initialized(self,base,expressions):
        # No automatic repeat of an expensive failed experiment.
        return self._worker_for(base).eval(expressions)

def snapshot():
    src=GP_SRC/"fatou_backend/vendor/fatou_fork.gp"
    engine=HERE/"engine_snapshot.gp"
    if not engine.exists():
        engine.write_bytes(src.read_bytes())
        (HERE/"provenance.json").write_text(json.dumps({
            "upstream":str(src),"sha256":hashlib.sha256(engine.read_bytes()).hexdigest(),
            "snapshot_date":"2026-10-12","upstream_and_shared_caches_modified":False},indent=2),encoding="utf-8")
    return engine

def engine(dps=30,extension=None,logged=False,timeout=120,even_extension=None):
    file=snapshot()
    if extension:
        content=snapshot().read_bytes()+b"\n"+Path(extension).read_bytes()
        even=Path(even_extension) if even_extension else HERE/"even8.gp"
        if even.exists():content+=b"\n"+even.read_bytes()
        file=HERE/("engine_with_hyper_"+hashlib.sha256(content).hexdigest()[:12]+".gp")
        if not file.exists():file.write_bytes(content)
    return (LoggedGP if logged else FatouGP)(gp_exe=find_default_gp_exe(),fatou_gp=file,dps=dps,
                   state_cache=False,init_timeout=150,eval_timeout=timeout,n_workers=1)

def stringify(v,digits=40):
    if isinstance(v,mp.mpc):
        return [mp.nstr(v.real,digits),mp.nstr(v.imag,digits)]
    return mp.nstr(v,digits)
