
CC = gcc
CXX = g++

CXXFLAGS = -Wall -Wno-deprecated

CFLAGS = -Wall			 	\
	-DG_DISABLE_DEPRECATED 	 	\
	-DGDK_DISABLE_DEPRECATED 	\
	-DGDK_PIXBUF_DISABLE_DEPRECATED \
	-DGTK_DISABLE_DEPRECATED

MODELS1 = model.h nonSpec.h nonSpec.cpp
#MODELS2 = model.h nonSpec.cpp
MODELS2 = nonSpec.cpp

BINARIES = handyDandy.o option2.o io_handler.o model.o nonSpec.o spatRes.o 
#BINARIES = option.o handyDandy.o io_handler.o model.o nonSpec.o 

all: ${BINARIES} applyModel

handyDandy.o: handyDandy.h handyDandy.cpp
	$(CXX) -c handyDandy.cpp $(CXXFLAGS)

option2.o: option2.h option2.cpp
	$(CXX) -c option2.cpp $(CXXFLAGS)

io_handler.o: io_handler.h io_handler.cpp handyDandy.o
	$(CXX) -c io_handler.cpp $(CXXFLAGS)

model.o: model.h model.cpp handyDandy.o
	$(CXX) -c model.cpp  $(CXXFLAGS)

nonSpec.o: nonSpec.h nonSpec.cpp model.o
	$(CXX) -c nonSpec.cpp $(CXXFLAGS)

spatRes.o: spatRes.h spatRes.cpp model.o
	$(CXX) -c spatRes.cpp $(CXXFLAGS)

applyModel: applyModel.cpp $(BINARIES)
	$(CXX) applyModel.cpp -o applyModel $(BINARIES) $(CXXFLAGS)

clean: 
	rm -f *.o *.gch applyModel
