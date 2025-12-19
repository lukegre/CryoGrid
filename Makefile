# --- Load Environment Variables ---
-include .env
export

# --- Configuration ---
# S3 paths
S3_PATH_PREFIX  := s3://spi-pamir-cryogrid/cryogrid_runs-luke
S3_PATH_FORCING := s3://spi-pamir-cryogrid/forcing/
# local paths on Euler
# if SCRATCH is not set, default to $(HOME)
WORKING_DIR     := $(or $(SCRATCH),$(HOME))
RUNS_DIR        := $(WORKING_DIR)/cryogrid-runs
FORCING_DIR     := $(HOME)/cryogrid-forcing/

# S3 related settings
export S3_ENDPOINT_URL := https://os.zhdk.cloud.switch.ch
export AWS_REQUEST_CHECKSUM_CALCULATION := WHEN_REQUIRED
export AWS_RESPONSE_CHECKSUM_VALIDATION := WHEN_REQUIRED

# Variable priority: 1. CLI (name=x) 2. .env (CRYOGRID_RUN_NAME)
RUN_NAME   := $(or $(name),$(CRYOGRID_RUN_NAME))
LOCAL_PATH := $(RUNS_DIR)/$(RUN_NAME)
S3_PATH    := $(S3_PATH_PREFIX)/$(RUN_NAME)/

.PHONY: help init install-aws download upload submit check-env check-name check-aws

help: ## Show this help message
	@echo "\033[1mUSAGE: make [target] [name=run-name]\033[0m"
	@grep -hE '^[a-zA-Z_-]+:.*?## .*$$' $(MAKEFILE_LIST) | sort | awk 'BEGIN {FS = ":.*?## "}; {printf "   \033[36m%-15s\033[0m  %s\n", $$1, $$2}'
	@echo "\033[1mCONFIGURATION\033[0m"
	@echo "   RUN_NAME         \033[1;33m$(RUN_NAME)\033[0m"
	@echo "   RUNS_DIR         \033[33m$(RUNS_DIR)\033[0m"
	@echo "   FORCING_DIR      \033[33m$(FORCING_DIR)\033[0m"
	@echo "   S3_PATH          \033[33m$(S3_PATH)\033[0m"
	@echo "   S3_PATH_FORCING  \033[33m$(S3_PATH_FORCING)\033[0m"
	

install-aws:  
	@if command -v aws >/dev/null 2>&1; then \
		echo "AWS CLI already installed. Skipping."; \
	else \
		echo "Downloading AWS CLI installer..."; \
		curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o awscliv2.zip; \
		unzip -q awscliv2.zip; \
		echo "Installing to $(HOME)/.local/bin..."; \
		./aws/install -i $(HOME)/.aws-cli -b $(HOME)/.local/bin --update; \
		echo "Installation complete. Ensure $(HOME)/.local/bin is in your PATH."; \
	fi
	
forcing: check-aws check-env check-name 
	@echo "Download ERA5 forcing for Pamirs to $(FORCING_DIR)..."
	@mkdir -p $(FORCING_DIR)
	@aws s3 sync $(S3_PATH_FORCING) $(FORCING_DIR) --endpoint-url $(S3_ENDPOINT_URL) --color on

dirs: 
	@echo "Creating necessary directories and adding symlink - $(HOME)/cryogrid-runs..."
	@mkdir -p "$(abspath $(RUNS_DIR))"
	@ln -snf "$(abspath $(RUNS_DIR))" "$(HOME)/cryogrid-runs"

init: dirs install-aws forcing  ## Setup scratch symlinks
	@echo "Initialization complete."
	
	
download-run: check-aws check-env check-name ## Sync files from S3 to local scratch
	@echo "Downloading $(S3_PATH)"
	@mkdir -p $(LOCAL_PATH)
	@aws s3 sync $(S3_PATH) $(LOCAL_PATH) --endpoint-url $(S3_ENDPOINT_URL) --color on

upload-run: check-aws check-env check-name ## Sync local scratch results to S3
	@echo "Uploading $(LOCAL_PATH)..."
	@aws s3 sync $(LOCAL_PATH) $(S3_PATH) --endpoint-url $(S3_ENDPOINT_URL) --color on

submit-run: check-name ## Submit the SLURM job
	@echo "Submitting $(RUN_NAME)..."
	@cd $(LOCAL_PATH) && sbatch sbatch_submit.sh

# --- Guards & Helpers ---
check-aws:
	@command -v aws >/dev/null 2>&1 || { echo >&2 "Error: AWS CLI not found. Run 'make install-aws' first."; exit 1; }

check-env:
ifndef AWS_ACCESS_KEY_ID
	$(error AWS_ACCESS_KEY_ID is missing in .env)
endif

check-name:
ifeq ($(RUN_NAME),)
	$(error RUN_NAME is not set. Use name=... or set CRYOGRID_RUN_NAME in .env)
endif