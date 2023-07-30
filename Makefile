COMPILER = g++
CFLAGS = -O3 -std=c++11 -I./minisat/

PRE_CFLAGS = ${CFLAGS} -c
TARGET = CAmpactor

SRC_DIR = src

BIN_DIR = bin

LSOPTIMIZER = calocalsearch
LSOPTIMIZER_TARGET = ${SRC_DIR}/${LSOPTIMIZER}.o
LSOPTIMIZER_CPP_FILE = ${SRC_DIR}/${LSOPTIMIZER}.cpp
LSOPTIMIZER_H_FILE = ${SRC_DIR}/${LSOPTIMIZER}.h
LSOPTIMIZER_SOURCE_FILES = ${LSOPTIMIZER_H_FILE} ${LSOPTIMIZER_CPP_FILE}

TARGET_FILES = ${LSOPTIMIZER_TARGET}
				
MAIN_SOURCE_FILE = ${SRC_DIR}/main.cpp

UPDATE = update
CLEAN = clean
CLEANUP = cleanup

all: ${TARGET_FILES} ${TARGET} ${UPDATE} ${CLEAN}

${LSOPTIMIZER_TARGET}: ${LSOPTIMIZER_SOURCE_FILES}
	${COMPILER} ${PRE_CFLAGS} ${LSOPTIMIZER_CPP_FILE} -o ${LSOPTIMIZER_TARGET}

${TARGET}: ${MAIN_SOURCE_FILE} ${TARGET_FILES}
	${COMPILER} ${CFLAGS} ${MAIN_SOURCE_FILE} ${TARGET_FILES} minisat/core/Solver.o -o ${TARGET}

${UPDATE}:
	chmod +x ${BIN_DIR}/*

${CLEAN}:
	rm -f *~
	rm -f ${SRC_DIR}/*.o
	rm -f ${SRC_DIR}/*~
	rm -f minisat/utils/*.or minisat/utils/*.o minisat/core/*.or minisat/core/*.o
	rm -f minisat/core/depend.mk

${CLEANUP}: ${CLEAN}
	rm -f ${TARGET}
