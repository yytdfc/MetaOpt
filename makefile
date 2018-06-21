
.PHONY: build clean update-include

build:
	mkdir -p build
	cd build; cmake ..
	make -C ./build -j8

test: build
	make test -C ./build

update-header:
	mkdir -p include
	wget https://raw.githubusercontent.com/arun11299/cpp-subprocess/master/subprocess.hpp -O ./include/subprocess.hpp
	wget https://raw.githubusercontent.com/mayah/tinytoml/master/include/toml/toml.h -O ./include/toml.hpp
	wget https://raw.githubusercontent.com/tanakh/cmdline/master/cmdline.h -O ./include/cmdline.hpp
	wget https://raw.githubusercontent.com/progschj/ThreadPool/master/ThreadPool.h -O ./include/thread_pool.hpp

clean:
	-rm -rf build
