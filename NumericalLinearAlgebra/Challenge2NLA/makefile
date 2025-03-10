INCLUDE_DIR = ./include
DATA_DIR = ./data
SRC_DIR = ./src

main: $(SRC_DIR)/* $(INCLUDE_DIR)/*
	g++ -c -I ${mkEigenInc} -o $(SRC_DIR)/utils.o $(SRC_DIR)/utils.cpp
	g++ -c -I ${mkEigenInc} -o $(SRC_DIR)/main.o $(SRC_DIR)/main.cpp
	g++ -o $(SRC_DIR)/main $(SRC_DIR)/utils.o $(SRC_DIR)/main.o

clean:
	rm -f $(SRC_DIR)/*.o
	rm -f $(SRC_DIR)/main

run:
	cd $(SRC_DIR) && ./main
