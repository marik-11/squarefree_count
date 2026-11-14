import subprocess
import time
from squarefree_count import squarefree_count, naive

print("n  naive_numpy  plain_numpy  plain_c  wheel_c")
for n in range(10, 19):
    x = 10 ** n
    a = "-"
    if n <= 15:
        t = time.time()
        v1 = naive(x)
        a = "%.3f" % (time.time() - t)
    t = time.time()
    v2 = squarefree_count(x)
    b = time.time() - t
    o1 = subprocess.run(["./sqf_plain", str(x)], capture_output=True, text=True).stdout.split()
    o2 = subprocess.run(["./sqf", str(x), "0.6"], capture_output=True, text=True).stdout.split()
    assert int(o1[0]) == v2 == int(o2[0])
    print(n, a, "%.3f" % b, o1[1], o2[1])
