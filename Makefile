all:
	$(MAKE) -C v21
	$(MAKE) -C v46
	$(MAKE) -C v70
	$(MAKE) -C v61
	$(MAKE) -C v60

clean:
	$(MAKE) -C v21 clean
	$(MAKE) -C v46 clean
	$(MAKE) -C v70 clean
	$(MAKE) -C v61 clean
	$(MAKE) -C v60 clean

.PHONY: all clean
