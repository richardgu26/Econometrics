MPICC     ?= mpic++
OFMPIINC  ?= $(shell $(MPICC) -showme:compile | sed -e "s/-pthread/-lpthread/g")
MPIINC    := $(OFMPIINC)
OFMPILIBS ?= $(shell $(MPICC) -showme:link | sed -e "s/-pthread/ /g")
MPILIBS   := $(OFMPILIBS)
MKOCTFILE ?= mkoctfile

all: MPI_Init.oct \
     MPI_Initialized.oct \
     MPI_Comm_rank.oct \
     MPI_Comm_size.oct \
     MPI_Finalize.oct \
     MPI_Finalized.oct \
     MPI_Send.oct \
     MPI_Recv.oct \
     MPI_Barrier.oct \
     MPI_Iprobe.oct \
     MPI_Probe.oct \
     MPI_Comm_Load.oct \
     MPI_Comm_Test.oct \
     MPI_Get_processor_name.oct \

%.oct: %.cc
	$(MKOCTFILE) $(MPIINC) $(MPILIBS) $<
clean:
	$(RM) *.o  *.oct *~

