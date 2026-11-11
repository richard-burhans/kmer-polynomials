"""Fetch the outputs of a laila_perf job whose driver died before fetching them (2026-11-11: the tunnel to laila's
Galaxy dropped mid-poll, so scripts/laila_perf.py exited without its fetch or its cleanup).

    python3 recon/kp10_vs_iblt_fetch.py <job id> <dest dir> [--cleanup]

Writes <dest>/run.log, counters.tsv, meta.tsv and job.json, the layout laila_perf.py writes, then prints the job's
state. With --cleanup, and only once the job is terminal, it purges the job's history (laila_perf's scratch history,
"laila-perf <label> (delete me)") and deletes the job's UDT, which the dead driver could not.
Uses cannabis-genome/scripts/common.py for the laila credentials.
"""

from __future__ import annotations

import json
import pathlib
import sys

sys.path.insert(0, "/home/rcb112/project/claude/bench/repos/cannabis-genome/scripts")
import common  # noqa: E402

TERMINAL = {"ok", "error", "deleted", "failed", "paused", "deleting"}
FETCH = {"counters": "counters.tsv", "log": "run.log", "meta": "meta.tsv"}


def main() -> int:
    import requests
    jid, dest = sys.argv[1], pathlib.Path(sys.argv[2])
    cleanup = "--cleanup" in sys.argv[3:]
    url, key = common.galaxy_credentials("laila")
    base, h = url.rstrip("/"), {"x-api-key": key}
    job = requests.get(f"{base}/api/jobs/{jid}", headers=h, params={"full": "true"}, timeout=120).json()
    state = job["state"]
    print(f"job {jid}: {state}")
    if state not in TERMINAL:
        return 2
    dest.mkdir(parents=True, exist_ok=True)
    rec = {"host": "laila", "job": jid, "state": state, "exit_code": job.get("exit_code"),
           "tool_id": job.get("tool_id"), "create_time": job.get("create_time"), "update_time": job.get("update_time"),
           "history_id": job.get("history_id"), "tool_stderr": (job.get("tool_stderr") or "")[-2000:],
           "fetched_by": "recon/kp10_vs_iblt_fetch.py (the laila_perf driver died before fetching)"}
    metrics = requests.get(f"{base}/api/jobs/{jid}/metrics", headers=h, timeout=120).json()
    if isinstance(metrics, list):
        rec["metrics"] = {m["name"]: m.get("raw_value", m.get("value")) for m in metrics
                          if not m["name"].startswith("processor_")}
    for name, fn in FETCH.items():
        ds = (job.get("outputs") or {}).get(name)
        if not ds:
            rec.setdefault("fetch_errors", {})[name] = "no such output"
            continue
        r = requests.get(f"{base}/api/datasets/{ds['id']}/display", headers=h, timeout=3600)
        if r.status_code != 200:
            rec.setdefault("fetch_errors", {})[name] = f"HTTP {r.status_code}"
            continue
        (dest / fn).write_bytes(r.content)
        print(f"  {fn}: {len(r.content):,} bytes")
    (dest / "job.json").write_text(json.dumps(rec, indent=1) + "\n")
    if cleanup:
        hist = job.get("history_id")
        if hist:
            name = requests.get(f"{base}/api/histories/{hist}", headers=h, timeout=60).json().get("name", "")
            if name.startswith("laila-perf ") and name.endswith("(delete me)"):
                r = requests.delete(f"{base}/api/histories/{hist}", headers=h, json={"purge": True}, timeout=120)
                print(f"  purged history {hist} ({name}): HTTP {r.status_code}")
            else:
                print(f"  history {hist} is not a laila-perf scratch history ({name!r}); left alone")
        tid = job.get("tool_uuid") or (job.get("tool_id") if str(job.get("tool_id", "")).count("-") == 4 else None)
        if tid:
            r = requests.delete(f"{base}/api/unprivileged_tools/{tid}", headers=h, timeout=60)
            print(f"  deleted UDT {tid}: HTTP {r.status_code}")
        else:
            print(f"  no UDT uuid on the job record (tool_id {job.get('tool_id')!r}); list it with GET /api/unprivileged_tools")
    return 0 if state == "ok" else 1


if __name__ == "__main__":
    sys.exit(main())
