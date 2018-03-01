EXE = qcf
OBJ = src/Bootstrap.o \
			src/MbRandom.o \
			src/SeqData.o \
			src/QCFTable.o \
			src/QCFData.o \
			src/Quartet.o \
			src/qcf.o

CXX = g++
CXXFLAGS = -Wall -O3 -g -std=c++11

.PHONY : clean test install uninstall

$(EXE) : $(OBJ)
	$(CXX) $(CXXFLAGS) -o $(EXE) $^
	@printf "\nTo install qcf system-wide run: sudo make install\n\n"

$(OBJ) : %.o: %.cpp
	$(CXX) $(CXXFLAGS) -o $@ -c $<

clean :
	@printf "\n Removing object (*.o) files and executable...\n\n"
	@rm -i $(OBJ) $(EXE)

test :
	@cd tests; qcf -i genes.txt -m maps.txt

install :
	@printf "\n Copying executable to /usr/local/bin...\n\n"
	@cp $(EXE) /usr/local/bin
	@printf " To uninstall, type: sudo make uninstall\n\n"

uninstall :
	@printf "\n Removing executable from /usr/local/bin...\n\n"
	@rm -i /usr/local/bin/$(EXE)
