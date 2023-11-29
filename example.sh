
# do the search changing the target one at a time 

# python search_sensor.py -d data/PSSM_logOddsPaper -k CHAK1 NEK11 YSK4 -t CDC7 -s 0.7 -o data/output/CDC7.txt -m 10
# python search_sensor.py -d data/PSSM_logOddsPaper -k CHAK1 NEK11 CDC7 -t YSK4 -s 0.7 -o data/output/YSK4.txt -m 10
# python search_sensor.py -d data/PSSM_logOddsPaper -k CHAK1 CDC7 YSK4 -t NEK11 -s 0.7 -o data/output/NEK11.txt -m 10
# python search_sensor.py -d data/PSSM_logOddsPaper -k CDC7 NEK11 YSK4 -t CHAK1 -s 0.7 -o data/output/CHAK1.txt -m 10


# # Define an array of file names
# files=("CDC7.txt" "YSK4.txt" "NEK11.txt" "CHAK1.txt")
# 
# 
# # Loop through each file in the array
# for file in "${files[@]}"
# do
#     concatenated_values=""
#     # Read the first column of each file, skip the first line, and concatenate the values
#     # Using awk: skip the first line (NR>1) and print the first column ($1)
#     for value in $(awk 'NR>1 {print $1}' data/output/"$file")
#     do
#         concatenated_values+="$value "
#     done
#     # Remove the trailing space
#     concatenated_values=${concatenated_values% }
#     echo $concatenated_values
#     #python plot_distribution.py -d data/PSSM_logOddsPaper -k CHAK1 NEK11 YSK4 -t CDC7 -s 0.7 -o data/output/CDC7.txt -m 10
#     #python plot_distribution.py -d data/PSSM_logOddsPaper -k CHAK1 NEK11 CDC7 -t YSK4 -s 0.7 -o data/output/YSK4.txt -m 10
#     #python plot_distribution.py -d data/PSSM_logOddsPaper -k CHAK1 CDC7 YSK4 -t NEK11 -s 0.7 -o data/output/NEK11.txt -m 10
#     #python plot_distribution.py -d data/PSSM_logOddsPaper -k CDC7 NEK11 YSK4 -t CHAK1 -s 0.7 -o data/output/CHAK1.txt -m 10
# done
# 
# # Remove the trailing space
# concatenated_values=${concatenated_values% }

# Pass the concatenated values to another command
# Replace 'your_command' with the actual command you want to use


# Define an array of elements
elements=("element1" "element2" "element3" "element4")
kinases=("CDC7" "YSK4" "NEK11" "CHAK1")

# Arrays to store commands
declare -a search_commands

# Prepare the arguments for each command
for (( i=0; i<${#kinases[@]}; i++ )); do
    target=${kinases[i]}
    background=()

    # Construct the background array
    for (( j=0; j<${#kinases[@]}; j++ )); do
        if [ $i -ne $j ]; then
            background+=("${kinases[j]}")
        fi
    done

    # Convert the background array to a space-separated string
    background_string="${background[*]}"

    # Store the command in the array
    search_commands+=("python search_sensor.py -d data/PSSM_logOddsPaper -k $background_string -t $target -s 0.7 -o data/output/'$target'.txt -m 10")
done

# Execute each search command in parallel
for cmd in "${search_commands[@]}"; do
    eval $cmd &  
done

wait

declare -a plot_commands
# Prepare the arguments for each command
for (( i=0; i<${#kinases[@]}; i++ )); do
    target=${kinases[i]}
    background=()

    # Construct the background array
    for (( j=0; j<${#kinases[@]}; j++ )); do
        if [ $i -ne $j ]; then
            background+=("${kinases[j]}")
        fi
    done

    # Convert the background array to a space-separated string
    background_string="${background[*]}"

    # Gather all the sequences
    concatenated_values=""
    # Read the first column of each file, skip the first line, and concatenate the values
    # Using awk: skip the first line (NR>1) and print the first column ($1)
    for value in $(awk 'NR>1 {print $1}' data/output/"$target".txt)
    do
        concatenated_values+="$value "
    done
    # Remove the trailing space
    concatenated_values=${concatenated_values% }

    # Store the command in the array
    plot_commands+=("python plot_distributions.py -d data/PSSM_logOddsPaper -k $background_string -t $target -s $concatenated_values -o data/output/'$target'.pdf")
done
# Wait for all background processes to finish

# Execute each search command in parallel
for cmd in "${plot_commands[@]}"; do
    eval $cmd &  
done

wait

echo "All commands have been executed."

