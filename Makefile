EXE = qcf
PY_EXE = scripts/qcf_boot.py
OBJ = src/Bootstrap.o \
			src/MbRandom.o \
			src/SeqData.o \
			src/QCFTable.o \
			src/QCFData.o \
			src/Quartet.o \
			src/qcf.o

CXX = g++
CXXFLAGS = -Wall -g -O3 -std=c++11
#CXXFLAGS += -Wno-sign-compare
#CXXFLAGS += -pedantic -Wextra

.PHONY : clean test install uninstall

$(EXE) : $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(EXE) $^
	@printf "\nTo install qcf system-wide run: sudo make install\n\n"

$(OBJ) : %.o: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

clean :
	@printf "\nRemoving object (*.o) files and executable...\n\n"
	@rm $(OBJ) $(EXE)

test :
	@printf "Running data in folder 'tests/'\n"
	@cd tests; ../qcf -i genes.txt -m map.txt --printRaw; ../scripts/qcf_boot.py -i out-raw.csv
	@printf "\n\nRunning data in folder 'example/'\n"
	@cd example; ../qcf -i genes.txt -m map.txt --printRaw; ../scripts/qcf_boot.py -i out-raw.csv

install :
	@printf "\nCopying executable to /usr/local/bin...\n\n"
	@cp $(EXE) $(PY_EXE) /usr/local/bin
	@printf "To uninstall, type: sudo make uninstall\n\n"

uninstall :
	@printf "\n Removing executable from /usr/local/bin...\n\n"
	@rm -i /usr/local/bin/$(EXE) /usr/local/bin/$(PY_EXE)
