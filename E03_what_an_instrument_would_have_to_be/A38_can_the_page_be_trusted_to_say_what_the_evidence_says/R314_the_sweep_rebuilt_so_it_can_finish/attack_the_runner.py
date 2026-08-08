r"""#875 · attack the RUNNER before trusting it — six vectors, each with the world it is aimed at

`#874` lost 109 minutes to a runner that was never attacked. `subprocess.run(capture_output=True,
timeout=...)` reads as safe and is not: **the timeout fires only when the pipe closes, and a
GRANDCHILD holds the pipe open after the child dies.** The failure is invisible from inside — the
parent sits at 0.0% CPU looking like slow work.

**P7: after building a lock, attack it — and after FIXING one, attack it again, because the fix is
where the new hole is.** Six vectors, each performed, each with its own expected verdict, and the
run FAILS LOUDLY if any of them does not behave. Every vector names the world it kills.

  V1 grandchild holds the output   -> the exact `#874` failure. Must be killed at the cap.
  V2 the grandchild must be DEAD   -> killing the child is not killing the group; check the pid.
  V3 a script that reads stdin     -> must not wait forever (`stdin=DEVNULL`).
  V4 a flood of stdout             -> must not deadlock (files, not pipes) and must be readable.
  V5 a clean exit                  -> the runner must not break the normal case (the negative side).
  V6 a fast failure                -> the classifier must still see stderr written to a file.

⚠ **V5 and V6 are the half that makes this a two-sided control.** A runner that killed everything
would pass V1–V4 and be useless; a control that only checks the dangerous cases is validated
against my imagination, not against the corpus.
"""
import os
import pathlib
import signal
import subprocess
import sys
import time

HERE = pathlib.Path(__file__).resolve().parent
STDIO = HERE / "results" / "_attack_stdio"
STDIO.mkdir(parents=True, exist_ok=True)
PY = sys.executable


def run_one(src, cap, tag):
    """The runner under attack — byte-identical in mechanism to the sweep's."""
    f = STDIO / f"{tag}.py"
    f.write_text(src)
    op, ep = STDIO / f"{tag}.out", STDIO / f"{tag}.err"
    t0 = time.time()
    with open(op, "wb") as so, open(ep, "wb") as se:
        p = subprocess.Popen([PY, str(f)], cwd=str(HERE), stdout=so, stderr=se,
                             stdin=subprocess.DEVNULL, start_new_session=True)
        timed_out = False
        while True:
            rc = p.poll()
            if rc is not None:
                break
            if time.time() - t0 > cap:
                timed_out = True
                try:
                    os.killpg(os.getpgid(p.pid), signal.SIGKILL)
                except (ProcessLookupError, PermissionError):
                    pass
                try:
                    rc = p.wait(timeout=30)
                except subprocess.TimeoutExpired:
                    rc = -9
                break
            time.sleep(0.2)
    return op.read_text(errors="replace"), ep.read_text(errors="replace"), rc, timed_out, time.time() - t0


CAP = 5
rows = []

# ---- V1/V2: a grandchild that outlives its parent and would hold a pipe forever
src = (f"import subprocess,sys,time,pathlib\n"
       f"g = subprocess.Popen([{PY!r}, '-c', 'import time;time.sleep(600)'])\n"
       f"pathlib.Path({str(STDIO / 'grandchild.pid')!r}).write_text(str(g.pid))\n"
       f"print('parent about to sleep', flush=True)\n"
       f"time.sleep(600)\n")
o, e, rc, to, dt = run_one(src, CAP, "v1")
gpid = int((STDIO / "grandchild.pid").read_text())
rows.append(("V1 grandchild holds the output — the exact `#874` failure",
             to and dt < CAP + 5, f"timed_out={to} after {dt:.1f}s (cap {CAP}s)"))
time.sleep(0.5)
alive = pathlib.Path(f"/proc/{gpid}").exists()
if alive:                                    # do not leave a 10-minute orphan behind either way
    try:
        os.kill(gpid, signal.SIGKILL)
    except ProcessLookupError:
        pass
rows.append(("V2 the GRANDCHILD is dead — killing the child is not killing the group",
             not alive, f"pid {gpid} alive after killpg: {alive}"))

# ---- V3: a script that waits on stdin
o, e, rc, to, dt = run_one("import sys\nprint('read:', sys.stdin.read())\n", CAP, "v3")
rows.append(("V3 a script that reads stdin must not wait forever (`stdin=DEVNULL`)",
             (not to) and rc == 0, f"rc={rc} timed_out={to} in {dt:.1f}s · stdout={o.strip()!r}"))

# ---- V4: a flood of stdout, which is what deadlocks a pipe
o, e, rc, to, dt = run_one("for i in range(200000): print('x'*80)\n", CAP, "v4")
rows.append(("V4 16 MB of stdout must not deadlock and must be readable",
             (not to) and rc == 0 and len(o) > 15_000_000,
             f"rc={rc} timed_out={to} bytes={len(o)} in {dt:.1f}s"))

# ---- V5: the normal case must still work (the other side of the control)
o, e, rc, to, dt = run_one("print('hello')\n", CAP, "v5")
rows.append(("V5 a clean script still exits 0 — a runner that kills everything passes V1-V4",
             rc == 0 and (not to) and o.strip() == "hello", f"rc={rc} stdout={o.strip()!r}"))

# ---- V6: stderr must survive to the file, or the classifier is blind
o, e, rc, to, dt = run_one("open('__no_such_file_v6__.csv')\n", CAP, "v6")
rows.append(("V6 stderr reaches the file — the classifier reads stderr, not the exit code alone",
             rc != 0 and "FileNotFoundError" in e, f"rc={rc} stderr_tail={e.strip()[-70:]!r}"))

print("=== ATTACK ON THE RUNNER — six vectors, each performed ===")
ok = True
for name, passed, detail in rows:
    ok &= bool(passed)
    print(f"  {'PASS' if passed else '**FAIL**'}  {name}\n         {detail}")
print(f"\n  => runner **{'SAFE TO USE' if ok else 'NOT SAFE — do not run the sweep'}**")
for p in STDIO.glob("*"):
    p.unlink()
STDIO.rmdir()
raise SystemExit(0 if ok else 2)
