#!/bin/bash

# Define the environment name and Python version
ENV_NAME="megs-v1"
PYTHON_VERSION="3.10" # You can change this to your desired Python version

# Function to create and activate environment using Conda
create_conda_environment() {
    conda create -n "$ENV_NAME" python="$PYTHON_VERSION" -y
    conda activate "$ENV_NAME"
}

# Function to create and activate environment using venv
create_venv_environment() {
    python3 -m venv "$ENV_NAME"
    source "$ENV_NAME/bin/activate"
}

# Check if Conda is available, if so use it, otherwise use venv
if command -v conda &>/dev/null; then
    echo "Conda is available. Creating a Conda environment..."
    create_conda_environment
else
    echo "Conda not available. Falling back to Python venv..."
    create_venv_environment
fi

echo "--------------------------------------------------"
# Install dependencies
echo "Installing dependencies from requirements.txt..."
python3 -m pip install -r requirements.txt

echo "--------------------------------------------------"
echo "Running setup.py install..."
python3 setup.py install --use

echo "--------------------------------------------------"
echo ""
echo "Setup complete. Virtual environment '$ENV_NAME' is ready."
echo "To activate this environment, use 'conda activate $ENV_NAME' (for Conda) or 'source $ENV_NAME/bin/activate' (for venv)."
echo "To deactivate this environment, use 'conda deactivate' (for Conda) or 'deactivate' (for venv)."
echo "--------------------------------------------------"