"""KP-13 Eisenstein revisit: pool kp13_fft_diag.py results over windows and compare encodings side by side.

Pooled as in results/ideas_spectral_2026-11-10.md: precision = hits / peaks over all windows; heavy-cluster recall =
found / heavy clusters; aligned-base-weighted recall = found weight / heavy weight; negative peaks per window.
    kp13_compare.py RESULT.json [RESULT.json ...]
"""
import json
import sys

for path in sys.argv[1:]:
    d = json.load(open(path))
    enc = d.get("encoding", "quarter")
    ws = d["windows"]
    for var in ws[0]["variants"]:
        for z in ("z6", "z10"):
            P = H = HC = F = W = FW = NP = 0
            for w in ws:
                v = w["variants"].get(var, {}).get(z)
                if v is None:
                    continue
                P += v["peaks"]
                H += (v["precision"] or 0) * v["peaks"]
                if w["heavy_clusters"]:
                    HC += w["heavy_clusters"]
                    F += (v["recall_heavy"] or 0) * w["heavy_clusters"]
                    W += w["heavy_weight"]
                    FW += (v["recall_heavy_weight"] or 0) * w["heavy_weight"]
                NP += v["neg_peaks"]
            print(f"{enc:10} {var:10} {z:4} peaks {P:4d}  precision {H / P if P else float('nan'):.3f}  "
                  f"heavy recall {F / HC:.3f}  weighted recall {FW / W:.3f}  neg/window {NP / len(ws):.2f}")
