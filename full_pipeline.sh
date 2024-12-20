#!/bin/bash

# This script runs the Python scripts plot_av.py and stats.py

check_environment() {
    if [[ "$(hostname)" == "MacBook-Pro-de-Nisso-6.local" ]]; then
        echo "Activating local Conda environment 'mconn'."
        conda init bash
        conda activate mconn
    else
        echo "Activating remote Conda environment 'mne'."
        conda init bash
        conda activate mne
    fi
}

# Activate the appropriate Conda environment
check_environment

# Run the Python scripts and the bash script
echo "Running plot_individuals.py.py..."
python3 plot_individuals.py
if [ $? -ne 0 ]; then
    echo "Error: plot_individuals.py failed."
    exit 1
fi

echo "Running plot_av.py..."
python3 plot_av.py
if [ $? -ne 0 ]; then
    echo "Error: plot_av.py failed."
    exit 1
fi

echo "Running stats.py..."
python3 stats.py
if [ $? -ne 0 ]; then
    echo "Error: stats.py failed."
    exit 1
fi

# if [[ "$(hostname)" == "MacBook-Pro-de-Nisso-6.local" ]]; then

#     echo "Running dl_results.sh..."
#     ./dl_results.sh
#     if [ $? -ne 0 ]; then
#         echo "Error: dl_results.sh failed."
#         exit 1
#     fi
# fi

# Structure of dl_results.sh:
# =================================
    # #!/bin/bash
    # Define the remote and local directories
    # REMOTE_USER=""        
    # REMOTE_HOST="" 
    # REMOTE_DIR="" 
    # LOCAL_DIR=""  
    # echo "Starting rsync transfer from $REMOTE_USER@$REMOTE_HOST:$REMOTE_DIR to $LOCAL_DIR"

    # # get the latest files
    # rsync -avz --progress "$REMOTE_USER@$REMOTE_HOST:$REMOTE_DIR" "$LOCAL_DIR"

    # # check if rsync was successful
    # if [ $? -eq 0 ]; then
    #     echo "Rsync completed successfully!"
    # else
    #     echo "Rsync encountered an error."
    # fi
# =================================