import random

def simulate_reads(n=10):
    bases = "ACGT"
    return ["".join(random.choice(bases) for _ in range(50)) for _ in range(n)]

reads = simulate_reads()

for r in reads:
    print(r)
