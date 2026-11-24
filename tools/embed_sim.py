#!/usr/bin/env python3
"""Replace the `const SIM={...};` table in index.html with the output of the Monte Carlo.

usage: python3 tools/embed_sim.py mc/eff.json [index.html]
"""
import json, re, sys

def fmt(v):
    return ("%.4g" % v) if v else "0"

def arr(rows):
    return "[" + ",".join("[" + ",".join(fmt(x) for x in r) + "]" for r in rows) + "]"

d = json.load(open(sys.argv[1]))
page = sys.argv[2] if len(sys.argv) > 2 else "index.html"
sim = "const SIM={E:%s,z:%s,N:%d" % (json.dumps(d["E"]).replace('"', ''), json.dumps(d["z"]), d["N"])
for k in ("nai", "bgo", "ge"):
    sim += ",%s:{pk:%s,tot:%s}" % (k, arr(d[k]["pk"]), arr(d[k]["tot"]))
sim += "};"
s = open(page).read()
s, n = re.subn(r"const SIM=\{.*?\}\};", lambda m: sim, s, count=1, flags=re.S)
assert n == 1, "SIM table not found"
open(page, "w").write(s)
print("updated", page, len(sim), "bytes")
