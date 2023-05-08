all:
	$(MAKE) -C v46
	$(MAKE) -C v21
	$(MAKE) -C vXXX

clean:
	$(MAKE) -C v46 clean
	$(MAKE) -C v21 clean
	$(MAKE) -C vXXX clean

.PHONY: all clean
