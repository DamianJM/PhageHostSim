# Phage host simulation environment
# Damian Magill - IFF 2026

import matplotlib
matplotlib.use("Agg")
import numpy as np
import pandas as pd
import random
import argparse
import matplotlib.pyplot as plt
#from matplotlib.animation import FuncAnimation, PillowWriter

# Simulation parameters  fed from bash file

def parse_args():
    parser = argparse.ArgumentParser(description="Bacterial growth and phage simulation")
    parser.add_argument('--width', type=int, default=100)
    parser.add_argument('--height', type=int, default=100)
    parser.add_argument('--initial_bacteria', type=int, default=50)
    parser.add_argument('--initial_phages', type=int, default=20)
    parser.add_argument('--carrying_capacity', type=int, default=500)
    parser.add_argument('--growth_rate', type=float, default=0.05)
    parser.add_argument('--diffusion_std', type=float, default=1.0)
    parser.add_argument('--time', type=int, default=50)
    parser.add_argument('--fps', type=int, default=5)
    parser.add_argument('--infection_radius', type=float, default=2.0)
    parser.add_argument('--burst_delay', type=int, default=10)
    parser.add_argument('--burst_size', type=int, default=5)
    parser.add_argument('--infection_probability', type=float, default=0.3)
    parser.add_argument('--chain_mode', type=str, default="True")
    parser.add_argument('--avg_chain_length', type=int, default=3)
    parser.add_argument('--max_chain_length', type=int, default=5)
    parser.add_argument('--output_file', type=str, default="simulation_output")
    parser.add_argument('--resistance_fraction', type=float, default=0.3)
    return parser.parse_args()

args = parse_args()

WIDTH, HEIGHT = args.width, args.height
INITIAL_BACTERIA, INITIAL_PHAGES = args.initial_bacteria, args.initial_phages
CARRYING_CAPACITY, GROWTH_RATE, DIFFUSION_STD = args.carrying_capacity, args.growth_rate, args.diffusion_std
TIME, FPS = args.time, args.fps
TOTAL_FRAMES = TIME * FPS
INFECTION_RADIUS, BURST_DELAY, BURST_SIZE = args.infection_radius, args.burst_delay, args.burst_size
INFECTION_PROBABILITY = args.infection_probability
CHAIN_MODE = args.chain_mode == "True"
AVG_CHAIN_LENGTH, MAX_CHAIN_LENGTH = args.avg_chain_length, args.max_chain_length
OUTPUT_FILE, RESISTANCE_FRACTION = args.output_file, args.resistance_fraction

# foundational functions for treating growth and interaction

def logistic_probability(N, r, K):
    return r * (1 - N / K)

def apply_repulsion(entities, threshold=2.0, strength=0.2):
    for i in range(len(entities)):
        for j in range(i + 1, len(entities)):
            dx, dy = entities[i].x - entities[j].x, entities[i].y - entities[j].y
            dist = np.hypot(dx, dy)
            if 0 < dist < threshold:
                repel = strength * (threshold - dist) / dist
                rx, ry = dx * repel, dy * repel
                entities[i].x += rx
                entities[i].y += ry
                entities[j].x -= rx
                entities[j].y -= ry
        if hasattr(entities[i], 'check_bounds'):
            entities[i].check_bounds()

# Phage and bacteria classes

class Bacterium:
    def __init__(self, x, y, sos):
        self.x, self.y = x, y
        self.infected = False
        self.infection_time = None
        self.sos = sos

    def check_bounds(self):
        self.x = max(0, min(self.x, WIDTH))
        self.y = max(0, min(self.y, HEIGHT))

    def move(self, dx=0, dy=0):
        self.x += dx
        self.y += dy
        self.check_bounds()

class BacteriaChain:
    def __init__(self, x, y, length, sos):
        self.bacteria = []
        self.angle = random.uniform(0, 2 * np.pi)
        for i in range(length):
            dx, dy = np.cos(self.angle) * i * 2, np.sin(self.angle) * i * 2
            self.bacteria.append(Bacterium(x + dx, y + dy, sos))

    def elongate(self, sos):
        """Extend the chain by adding a daughter cell at the growing tip."""
        tip = self.bacteria[-1]
        nx = tip.x + np.cos(self.angle) * 2
        ny = tip.y + np.sin(self.angle) * 2
        new_bact = Bacterium(nx, ny, sos)
        new_bact.x = max(0, min(nx, WIDTH))
        new_bact.y = max(0, min(ny, HEIGHT))
        self.bacteria.append(new_bact)

    def split(self):
        """Split chain at the midpoint; return the new daughter chain."""
        mid = len(self.bacteria) // 2
        daughter = BacteriaChain.__new__(BacteriaChain)
        daughter.bacteria = self.bacteria[mid:]
        # Slight angle divergence between parent and daughter chain
        daughter.angle = (self.angle + random.uniform(-0.3, 0.3)) % (2 * np.pi)
        self.angle = (self.angle + random.uniform(-0.3, 0.3)) % (2 * np.pi)
        self.bacteria = self.bacteria[:mid]
        return daughter

    def move(self):
        dx, dy = np.random.normal(0, DIFFUSION_STD, 2)
        for b in self.bacteria:
            b.move(dx, dy)

    def center(self):
        xs, ys = [b.x for b in self.bacteria], [b.y for b in self.bacteria]
        return np.mean(xs), np.mean(ys)

class Phage:
    def __init__(self, x, y, sos):
        self.x, self.y, self.viable, self.sos = x, y, True, sos

    def check_bounds(self):
        self.x = max(0, min(self.x, WIDTH))
        self.y = max(0, min(self.y, HEIGHT))

    def move(self):
        dx, dy = np.random.normal(0, DIFFUSION_STD, 2)
        self.x += dx
        self.y += dy
        self.check_bounds()

    def change_viability(self):
        self.viable = not self.viable

# Initialize bacteria chains and phages
margin = MAX_CHAIN_LENGTH * 2
bacteria = []
while sum(len(c.bacteria) for c in bacteria) < INITIAL_BACTERIA:
    length = min(max(1, int(np.random.exponential(AVG_CHAIN_LENGTH))), MAX_CHAIN_LENGTH)
    sos = 0 if random.random() < RESISTANCE_FRACTION else 1
    x, y = random.uniform(margin, WIDTH-margin), random.uniform(margin, HEIGHT-margin)
    bacteria.append(BacteriaChain(x, y, length, sos))

phages = [Phage(random.uniform(0, WIDTH), random.uniform(0, HEIGHT), 1) for _ in range(INITIAL_PHAGES)]
population_over_time = []

fig, ax = plt.subplots(figsize=(6, 6))
ax.set_xlim(0, WIDTH)
ax.set_ylim(0, HEIGHT)
sc_bacteria = ax.scatter([], [], c=[], cmap='coolwarm', s=20)
sc_phages = ax.scatter([], [], c='black', s=5)

def update(frame):
    global bacteria, phages
    new_chains, new_phages = [], []
    total_bacteria = sum(len(c.bacteria) for c in bacteria)
    contact_count, contact_sensitive, contact_resistant, bb_contact_count = 0, 0, 0, 0

    for chain in bacteria:
        chain.move()
        for b in chain.bacteria[:]:
            if b.infected and frame - b.infection_time >= BURST_DELAY:
                new_phages.extend([Phage(b.x, b.y, b.sos) for _ in range(BURST_SIZE)])
                chain.bacteria.remove(b)

        if total_bacteria < CARRYING_CAPACITY:
            if chain.bacteria and random.random() < logistic_probability(total_bacteria, GROWTH_RATE, CARRYING_CAPACITY):
                sos = chain.bacteria[0].sos
                chain.elongate(sos)
                total_bacteria += 1
                # Split when the chain exceeds max length (binary fission)
                if len(chain.bacteria) > MAX_CHAIN_LENGTH:
                    new_chains.append(chain.split())

    bacteria = [c for c in bacteria if c.bacteria]
    bacteria.extend(new_chains)
    all_bacteria = [b for chain in bacteria for b in chain.bacteria]

    for i, b1 in enumerate(all_bacteria):
        for b2 in all_bacteria[i+1:]:
            if np.hypot(b1.x - b2.x, b1.y - b2.y) < INFECTION_RADIUS:
                bb_contact_count += 1

    apply_repulsion(all_bacteria, threshold=INFECTION_RADIUS, strength=0.2)

    for p in phages:
        p.move()
        for b in all_bacteria:
            if np.hypot(b.x - p.x, b.y - p.y) < INFECTION_RADIUS:
                contact_count += 1
                if b.sos == 1:
                    contact_sensitive += 1
                else:
                    contact_resistant += 1
                if not b.infected and random.random() < INFECTION_PROBABILITY and p.sos == b.sos:
                    b.infected, b.infection_time = True, frame
                    p.change_viability()
                    break

    apply_repulsion(phages, threshold=INFECTION_RADIUS, strength=0.2)
    phages = [p for p in phages if p.viable]
    phages.extend(new_phages)

    if all_bacteria:
        sc_bacteria.set_offsets([[b.x, b.y] for b in all_bacteria])
        sc_bacteria.set_array(np.array([b.sos for b in all_bacteria]))
    else:
        sc_bacteria.set_offsets(np.empty((0, 2)))
        sc_bacteria.set_array(np.array([]))

    if phages:
        sc_phages.set_offsets([[p.x, p.y] for p in phages])
    else:
        sc_phages.set_offsets(np.empty((0, 2)))

    infected_count = sum(b.infected for b in all_bacteria)
    total_sensitive = sum(b.sos == 1 for b in all_bacteria)
    total_resistant = sum(b.sos == 0 for b in all_bacteria)

    population_over_time.append({
        'Frame': frame,
        'Total_Bacteria': len(all_bacteria),
        'Infected_Bacteria': infected_count,
        'Sensitive_Bacteria': total_sensitive,
        'Resistant_Bacteria': total_resistant,
        'Phages': len(phages),
        'Phage_Bacteria_Contacts': contact_count,
        'Phage_Sensitive_Contacts': contact_sensitive,
        'Phage_Resistant_Contacts': contact_resistant,
        'Bacteria_Bacteria_Contacts': bb_contact_count
    })

for frame in range(TOTAL_FRAMES):
    update(frame)

#ani = FuncAnimation(fig, update, frames=TOTAL_FRAMES, interval=1000/FPS, repeat=False)
#ani.save(f"{OUTPUT_FILE}.gif", writer=PillowWriter(fps=FPS))
df = pd.DataFrame(population_over_time)
df.to_csv(f'{OUTPUT_FILE}.csv', index=False)
print(f"Simulation complete — GIF saved as {OUTPUT_FILE}.gif")

