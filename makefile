pr.exe : main.cpp
	@g++ -std=c++11 -O3 $^ -o $@

run:
	@./pr.exe $(filter-out $@,$(MAKECMDGOALS))

%:
	@:

compile : pr.exe
.PHONY : compile
.SILENT :