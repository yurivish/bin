.PHONY: wasm
wasm:
	zig cc -o ../dist/bitvector.wasm bitvector.c walloc.c -shared -target wasm32-freestanding \
	 -Wall -Wextra -Wpedantic -Werror \
	 -std=c17 -nostdlib \
	 -g0 \
	 -O3

# note: -Oz reduces code size 
# useful tool: https://evmar.github.io/weave/
