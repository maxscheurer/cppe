CURRENT_DIR=$(PWD)
DALTON_FULL="dalton_full"
DALTON_IF="dalton_gen1int_interface"
# Checkout or update dalton_full folder
if [ -d "$DALTON_FULL" ]; then
	cd $DALTON_FULL
	git fetch
	git pull
else
	git clone https://gitlab.com/dalton/dalton.git dalton_full
fi
cd $CURRENT_DIR
# Copy needed files from dalton_full to dalton_gen1int_interface

if [ -d "$DALTON_IF" ]; then
	cp -r $DALTON_FULL/DALTON/gen1int $DALTON_IF/
else
	echo "cppe repository not checked out correctly"
fi
