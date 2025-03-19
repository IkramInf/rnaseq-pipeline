#!/bin/bash

# RNA-seq Pipeline Installation Script
# This script installs all necessary dependencies for the RNA-seq pipeline
# Version: 1.0

set -e  # Exit on error

# Check for sudo access
check_sudo() {
    if ! sudo -v &> /dev/null; then
        echo "Error: This script requires sudo privileges to install dependencies."
        echo "Please run with sudo or as root."
        exit 1
    fi
}

# Detect OS type
detect_os() {
    if [ -f /etc/os-release ]; then
        . /etc/os-release
        OS=$NAME
        
        if [[ "$OS" == *"Ubuntu"* ]] || [[ "$OS" == *"Debian"* ]]; then
            PACKAGE_MANAGER="apt-get"
        elif [[ "$OS" == *"CentOS"* ]] || [[ "$OS" == *"Red Hat"* ]] || [[ "$OS" == *"Fedora"* ]]; then
            PACKAGE_MANAGER="yum"
        else
            echo "Warning: Unsupported OS. This script is designed for Ubuntu/Debian or CentOS/RHEL/Fedora."
            echo "Attempting to continue with apt-get..."
            PACKAGE_MANAGER="apt-get"
        fi
    else
        echo "Warning: Unable to determine OS. This script is designed for Ubuntu/Debian or CentOS/RHEL/Fedora."
        echo "Attempting to continue with apt-get..."
        PACKAGE_MANAGER="apt-get"
    fi
    
    echo "Detected OS: $OS"
    echo "Using package manager: $PACKAGE_MANAGER"
}

# Check for existing tools
check_installed() {
    local cmd="$1"
    if command -v "$cmd" &> /dev/null; then
        echo "✓ $cmd is already installed"
        return 0
    else
        echo "✗ $cmd is not installed"
        return 1
    fi
}

# Install system dependencies
install_system_dependencies() {
    echo "Updating package lists..."
    if [ "$PACKAGE_MANAGER" = "apt-get" ]; then
        sudo apt-get update -qq
        
        echo "Installing system dependencies..."
        sudo apt-get install -y \
            wget \
            curl \
            default-jre \
            build-essential \
            zlib1g-dev \
            libncurses5-dev \
            libbz2-dev \
            liblzma-dev \
            libcurl4-openssl-dev \
            libssl-dev \
            python3 \
            python3-pip \
            python3-venv
    else
        sudo yum update -y -q
        
        echo "Installing system dependencies..."
        sudo yum install -y \
            wget \
            curl \
            java-11-openjdk \
            make \
            gcc \
            gcc-c++ \
            zlib-devel \
            ncurses-devel \
            bzip2-devel \
            xz-devel \
            libcurl-devel \
            openssl-devel \
            python3 \
            python3-pip \
            python3-devel
    fi
    
    echo "System dependencies installed"
}

# Install Python packages
install_python_packages() {
    echo "Installing Python packages..."
    
    # Create virtual environment
    python3 -m venv ~/rna_seq_env
    source ~/rna_seq_env/bin/activate
    
    # Install packages
    pip install --upgrade pip
    pip install multiqc fastqc-tools cutadapt biopython pandas matplotlib seaborn
    
    echo "Python packages installed in virtual environment: ~/rna_seq_env"
    echo "To activate: source ~/rna_seq_env/bin/activate"
}

# Install STAR aligner
install_star() {
    if check_installed "STAR"; then
        return
    fi
    
    echo "Installing STAR aligner..."
    
    # Create temporary directory
    mkdir -p ~/temp_install
    cd ~/temp_install
    
    # Download and install STAR
    wget https://github.com/alexdobin/STAR/archive/2.7.10b.tar.gz
    tar -xzf 2.7.10b.tar.gz
    cd STAR-2.7.10b/source
    make -j$(nproc)
    
    # Move to /usr/local/bin
    sudo mv STAR /usr/local/bin/
    
    # Clean up
    cd ~
    rm -rf ~/temp_install
    
    echo "STAR aligner installed"
}

# Install Trimmomatic
install_trimmomatic() {
    if check_installed "trimmomatic"; then
        return
    fi
    
    echo "Installing Trimmomatic..."
    
    # Create directory
    mkdir -p ~/tools/trimmomatic
    cd ~/tools/trimmomatic
    
    # Download Trimmomatic
    wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.39.zip
    unzip Trimmomatic-0.39.zip
    
    # Create wrapper script
    cat > ~/tools/trimmomatic/trimmomatic << EOF
#!/bin/bash
java -jar ~/tools/trimmomatic/Trimmomatic-0.39/trimmomatic-0.39.jar "\$@"
EOF
    
    chmod +x ~/tools/trimmomatic/trimmomatic
    sudo ln -sf ~/tools/trimmomatic/trimmomatic /usr/local/bin/
    
    echo "Trimmomatic installed"
}

# Check installation
verify_installation() {
    echo "Verifying installation..."
    
    tools=("STAR" "trimmomatic" "fastqc" "multiqc" "java")
    all_installed=true
    
    for tool in "${tools[@]}"; do
        if ! check_installed "$tool"; then
            all_installed=false
        fi
    done
    
    if [ "$all_installed" = true ]; then
        echo "All tools are installed correctly!"
    else
        echo "Warning: Some tools might not be installed properly. Please check the error messages above."
    fi
}

# Create activation script
create_activation_script() {
    cat > ~/activate_rna_seq_pipeline.sh << EOF
#!/bin/bash
# RNA-seq Pipeline Environment Activation Script

# Activate Python virtual environment
source ~/rna_seq_env/bin/activate

# Set PATH
export PATH=\$PATH:~/tools/trimmomatic

# Print activation message
echo "RNA-seq pipeline environment activated!"
echo "Run ./rna_seq_pipeline.sh --help for usage information"
EOF

    chmod +x ~/activate_rna_seq_pipeline.sh
    echo "Created activation script: ~/activate_rna_seq_pipeline.sh"
    echo "Source this script before running the pipeline: source ~/activate_rna_seq_pipeline.sh"
}

# Main function
main() {
    echo "==== RNA-seq Pipeline Installation ===="
    
    # Check for sudo access
    check_sudo
    
    # Detect OS
    detect_os
    
    # Install system dependencies
    install_system_dependencies
    
    # Install Python packages
    install_python_packages
    
    # Install tools
    install_star
    install_trimmomatic
    
    # Verify installation
    verify_installation
    
    # Create activation script
    create_activation_script
    
    echo "==== Installation Complete ===="
    echo
    echo "To use the RNA-seq pipeline:"
    echo "1. Source the activation script: source ~/activate_rna_seq_pipeline.sh"
    echo "2. Run the pipeline: ./rna_seq_pipeline.sh --help"
    echo
    echo "Thank you for installing the RNA-seq pipeline!"
}

# Execute main function
main
