import yaml

# Paths to the reformatted files
reformatted_r1 = "dbit/reformatted.R1.fastq.gz"
reformatted_r2 = "dbit/reformatted.R2.fastq.gz"

# Load the YAML configuration file
config_path = "config.yaml"
with open(config_path, 'r') as file:
    config = yaml.safe_load(file)

# Update the paths for r1 and r2
config['project']['r1'] = reformatted_r1
config['project']['r2'] = reformatted_r2

# Save the updated configuration back to the YAML file
with open(config_path, 'w') as file:
    yaml.dump(config, file)

print(f"Updated config['project']['r1'] to {reformatted_r1}")
print(f"Updated config['project']['r2'] to {reformatted_r2}")
