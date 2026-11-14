import random
import subprocess
import sys
from squarefree_count import squarefree_count, naive

known = [1, 7, 61, 608, 6083, 60794, 607926, 6079291, 60792694, 607927124, 6079270942,
         60792710280, 607927102274, 6079271018294, 60792710185947, 607927101854103,
         6079271018540405, 60792710185403794, 607927101854022750, 6079271018540280875,
         60792710185402613302, 607927101854026645617]
bad = 0
for n in range(0, 17):
    v = squarefree_count(10 ** n) if n >= 2 else naive(10 ** n)
    if v != known[n]:
        bad += 1
        print("py mismatch", n, v, known[n])
for prog in ("./sqf_plain", "./sqf"):
    for n in range(0, 22):
        out = subprocess.run([prog, str(10 ** n)], capture_output=True, text=True).stdout.split()[0]
        if int(out) != known[n]:
            bad += 1
            print(prog, "mismatch", n, out, known[n])
for _ in range(300):
    x = random.randint(1, 10 ** 12)
    c = random.choice([0.05, 0.3, 1.0, 3.0])
    a = subprocess.run(["./sqf_plain", str(x), str(c)], capture_output=True, text=True).stdout.split()[0]
    b = subprocess.run(["./sqf", str(x), str(c)], capture_output=True, text=True).stdout.split()[0]
    if a != b:
        bad += 1
        print("wheel mismatch", x, c, a, b)
for _ in range(300):
    x = random.randint(1, 10 ** 7)
    D = random.randint(2, 4000)
    if squarefree_count(x, D) != naive(x):
        bad += 1
        print("random mismatch", x, D)
print("ok" if bad == 0 else "fail", bad)
