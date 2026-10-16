#!/bin/bash

# Default parameters
WIDTH=100
HEIGHT=100
INITIAL_BACTERIA=30
CARRYING_CAPACITY=200
GROWTH_RATE=0.05
DIFFUSION_STD=1.0
TIME=20
FPS=20
INFECTION_RADIUS=2
CHAIN_MODE=True

# Parameter sweeps
INITIAL_PHAGES_LIST=(1 5 10 20 30 60)
MAX_CHAIN_LENGTH_LIST=(4 8 14 20)
INFECTION_PROBABILITY_LIST=(0.7)

BURST_SIZE_LIST=(10)
BURST_DELAY_LIST=(20)
RESISTANCE_FRACTION_LIST=(0)

NUMBER_REPLICATES=30

# Parallelisation
MAX_JOBS=1
job_count=0

# Compute total number of simulations 

TOTAL_SIMS=$(( \
${#INITIAL_PHAGES_LIST[@]} * \
${#MAX_CHAIN_LENGTH_LIST[@]} * \
${#INFECTION_PROBABILITY_LIST[@]} * \
${#BURST_SIZE_LIST[@]} * \
${#BURST_DELAY_LIST[@]} * \
${#RESISTANCE_FRACTION_LIST[@]} * \
NUMBER_REPLICATES ))

sim_counter=0

run_simulation () {
    phage_count=$1
    chain_length=$2
    infection_prob=$3
    burst_size=$4
    burst_delay=$5
    resistance_fraction=$6
    rep=$7
    sim_id=$8

    # Clean float values for filenames
    ip_tag=$(echo "$infection_prob" | sed 's/\./p/g')
    rf_tag=$(echo "$resistance_fraction" | sed 's/\./p/g')

    output_name="sim_p${phage_count}_c${chain_length}_ip${ip_tag}_bs${burst_size}_bd${burst_delay}_rf${rf_tag}_r${rep}"

    echo "[${sim_id}/${TOTAL_SIMS}] Running: phage=$phage_count chain=$chain_length ip=$infection_prob bs=$burst_size bd=$burst_delay rf=$resistance_fraction rep=$rep"

    python3 chain_length_phage_simulation.py \
        --width "$WIDTH" \
        --height "$HEIGHT" \
        --initial_bacteria "$INITIAL_BACTERIA" \
        --initial_phages "$phage_count" \
        --carrying_capacity "$CARRYING_CAPACITY" \
        --growth_rate "$GROWTH_RATE" \
        --diffusion_std "$DIFFUSION_STD" \
        --time "$TIME" \
        --fps "$FPS" \
        --infection_radius "$INFECTION_RADIUS" \
        --burst_delay "$burst_delay" \
        --burst_size "$burst_size" \
        --infection_probability "$infection_prob" \
        --chain_mode "$CHAIN_MODE" \
        --avg_chain_length "$chain_length" \
        --max_chain_length "$chain_length" \
        --output_file "$output_name" \
        --resistance_fraction "$resistance_fraction" \
        > "${output_name}.log" 2>&1
}

# Main loop
for phage_count in "${INITIAL_PHAGES_LIST[@]}"; do
    for chain_length in "${MAX_CHAIN_LENGTH_LIST[@]}"; do
        for infection_prob in "${INFECTION_PROBABILITY_LIST[@]}"; do
            for burst_size in "${BURST_SIZE_LIST[@]}"; do
                for burst_delay in "${BURST_DELAY_LIST[@]}"; do
                    for resistance_fraction in "${RESISTANCE_FRACTION_LIST[@]}"; do
                        for rep in $(seq 1 $NUMBER_REPLICATES); do

                            ((sim_counter++))

                            run_simulation \
                                "$phage_count" \
                                "$chain_length" \
                                "$infection_prob" \
                                "$burst_size" \
                                "$burst_delay" \
                                "$resistance_fraction" \
                                "$rep" \
                                "$sim_counter" &

                            ((job_count++))

                            # Control parallel jobs
                            if [ "$job_count" -ge "$MAX_JOBS" ]; then
                                wait
                                job_count=0
                            fi

                        done
                    done
                done
            done
        done
    done
done

wait
echo "All simulations complete: ${TOTAL_SIMS} runs"