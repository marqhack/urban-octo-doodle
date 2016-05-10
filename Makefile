CC = g++ -g -std=c++11 -pedantic

CC_DEBUG = @$(CC) -pedantic
CC_RELEASE = @$(CC) -O3 -DNDEBUG

G_SRC = src/*.cpp *.cpp

# need libpng to build
#
G_INC = -Iinclude -Iapps -I/opt/local/include -L/opt/local/lib

all: image tests

image : $(G_SRC) apps/image.cpp apps/image_recs.cpp
	$(CC_DEBUG) $(G_INC) $(G_SRC) apps/image.cpp apps/image_recs.cpp -lpng -o image

tests : $(G_SRC)  apps/tests.cpp apps/test_recs.cpp
	$(CC_DEBUG) $(G_INC) $(G_SRC) apps/tests.cpp apps/test_recs.cpp -lpng -o tests

# needs xwindows to build
#
X_INC = -I/opt/X11/include -L/opt/X11/lib

DRAW_SRC = apps/draw.cpp apps/GWindow.cpp
draw: $(DRAW_SRC) $(G_SRC)
	$(CC_RELEASE) $(X_INC) $(G_INC) $(G_SRC) $(DRAW_SRC) -lpng -lX11 -o draw


clean:
	@rm -rf image tests draw *.png *.dSYM

