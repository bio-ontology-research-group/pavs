# Create a directory for the tool
mkdir -p ~/tools/phenopacket-tools
cd ~/tools/phenopacket-tools

# Download the latest release
wget https://github.com/phenopackets/phenopacket-tools/releases/latest/download/phenopacket-tools-cli.jar

# Make it executable (optional, for convenience)
echo 'alias phenopacket-tools="java -jar ~/tools/phenopacket-tools/phenopacket-tools-cli.jar"' >> ~/.bashrc
source ~/.bashrc

# Test it
phenopacket-tools --help
