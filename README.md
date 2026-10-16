# PhageHostSim

A simple spatial simulation of bacterial growth and bacteriophage infection. The model represents bacteria as chains, allows sensitive and resistant cells, and records population and contact statistics over time.

## Files

- `chain_length_phage_simulation.py` – runs a single simulation.
- `run_simulation.sh` – runs parameter sweeps with multiple replicates.

## Requirements

- Python 3
- NumPy
- pandas
- Matplotlib
- Bash, for running the parameter sweep

Install the Python dependencies with:

```bash
pip install numpy pandas matplotlib
```

## Usage

### Run one simulation

```bash
python3 chain_length_phage_simulation.py --output_file simulation_output
```

Use `--help` to view all available parameters:

```bash
python3 chain_length_phage_simulation.py --help
```

### Run the parameter sweep

Make the shell script executable, then run it from the same directory as the Python script:

```bash
chmod +x run_simulation.sh
./run_simulation.sh
```

The sweep settings, replicates, and maximum number of parallel jobs can be edited near the top of `run_simulation.sh`.

## Output

Each simulation produces:

- A CSV file containing bacterial populations, infected-cell counts, phage counts, and contact statistics for every frame.
- A log file when the simulation is launched through `run_simulation.sh`.

Output filenames encode the main parameter values and replicate number.

## Notes

The simulation is stochastic, so results can differ between replicates. GIF generation is currently disabled in the Python script.
