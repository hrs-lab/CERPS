"""Run all simulation batches (checkpointed, 2 workers).

Usage:  python run_batches.py [name_filter ...] [--force]
"""
import os
import sys
import time
import traceback

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)

from common import build_jobs, job_path, DATA  # noqa: E402


def run_one(name, fn):
    t0 = time.time()
    out = fn()
    out["__name__"] = name
    tmp = job_path(name) + ".tmp.npz"
    np.savez_compressed(tmp, **out)
    os.replace(tmp, job_path(name))
    return time.time() - t0


def main():
    argv = list(sys.argv[1:])
    budget = 400.0
    if "--budget" in argv:
        k = argv.index("--budget")
        budget = float(argv[k + 1])
        del argv[k:k + 2]
    args = [a for a in argv if not a.startswith("--")]
    force = "--force" in sys.argv
    jobs = build_jobs()
    todo = []
    for name in jobs:
        if args and not any(a in name for a in args):
            continue
        if not force and os.path.exists(job_path(name)):
            continue
        todo.append(name)
    print(f"{ len(jobs) } jobs defined, { len(todo) } to run", flush=True)

    order = {n: i for i, n in enumerate([
        "base", "fig8", "fig9", "fig11", "fig4", "fig12", "fig7",
        "fig3", "fig10", "fig11c"])}
    todo.sort(key=lambda n: next(
        (v for k, v in order.items() if n.startswith(k)), 99))

    import time as _t
    from concurrent.futures import ProcessPoolExecutor, as_completed
    log = open(os.path.join(DATA, "_runlog.txt"), "a")
    t_start = _t.time()
    stop = False
    with ProcessPoolExecutor(max_workers=2) as ex:
        futs = {}
        for n in todo:
            if _t.time() - t_start > 40:   # only start jobs safe to finish
                break
            futs[ex.submit(run_one, n, jobs[n])] = n
        for fut in as_completed(futs):
            name = futs[fut]
            try:
                dt = fut.result()
                msg = f"{ name:28s} done in { dt:7.1f}s"
            except Exception:
                msg = f"{ name:28s} FAILED\n{ traceback.format_exc() }"
            print(msg, flush=True)
            log.write(msg + "\n")
            log.flush()
            if _t.time() - t_start > budget:
                print("budget reached, stopping submissions", flush=True)
                break
    log.close()
    print(f"chunk finished: {len(futs)} jobs submitted, "
          f"{_t.time() - t_start:.0f}s elapsed", flush=True)


if __name__ == "__main__":
    main()
