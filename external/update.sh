CURRENT_DIR=$(PWD)
DALTON_FULL="dalton_full"
# Checkout or update dalton_full folder
if [ -d "$DALTON_FULL" ]; then
	cd $DALTON_FULL
	git fetch
	git pull
else
	git clone https://gitlab.com/dalton/dalton.git dalton_full
fi

# Copy needed files from dalton_full to dalton_gen1int_interface
