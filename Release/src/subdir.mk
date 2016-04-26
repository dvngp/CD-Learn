################################################################################
# Automatically-generated file. Do not edit!
################################################################################

# Add inputs and outputs from these tool invocations to the build variables 
CPP_SRCS += \
../src/BE.cpp \
../src/Function.cpp \
../src/GM.cpp \
../src/MAIN.cpp \
../src/Variable.cpp

OBJS += \
./src/BE.o \
./src/Function.o \
./src/GM.o \
./src/MAIN.o \
./src/Variable.o 

CPP_DEPS += \
./src/BE.d \
./src/Function.d \
./src/GM.d \
./src/MAIN.d \
./src/Variable.d 

#CFLAGS = -D _LOGGING_ -D _SOFTEVID_
#CFLAGS = -D _LOGGING_
CFLAGS = -D _EXACTCOUNT_
# Each subdirectory must supply rules for building sources it contributes
src/%.o: ../src/%.cpp
	@echo 'Building file: $<'
	@echo 'Invoking: GCC C++ Compiler'
	g++ $(CFLAGS) -static -DNDEBUG -O3 -c -fmessage-length=0 -MMD -MP -MF"$(@:%.o=%.d)" -MT"$(@:%.o=%.d)" -o "$@" "$<"
	@echo 'Finished building: $<'
	@echo ' '


