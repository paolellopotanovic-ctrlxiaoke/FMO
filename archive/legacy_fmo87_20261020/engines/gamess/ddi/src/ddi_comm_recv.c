/* -------------------------------------------------------------------- *\
 * Distributed Data Interface
 * ==========================
 * 
 * Global objects
 * 
 * Author: Ryan M. Olson
 * 18 Aug 10 - MWS - correct arg checking messages, to not say 'send'!
\* -------------------------------------------------------------------- */
 # include "ddi_base.h"


/* -------------------------------------------------------------------- *\
   DDI_Recv(bufer,size,from)
   =========================
   [OUT] buff - address in which incoming data is placed.
   [IN]  size - size of incoming data stream.
   [IN]  from - rank of process sending the data stream.
   
   Used for point-to-point communications and large internal data
   tranfers to/from data servers.  Using working communicator.
\* -------------------------------------------------------------------- */
   void DDI_Recv(void *buff,size_t size,int from) {
      const DDI_Comm *comm = (DDI_Comm *) Comm_find(DDI_WORKING_COMM);
      DEBUG_OUT(LVL1,(stdout,"%s: entering DDI_Recv.\n",DDI_Id()))
      Comm_recv(buff,size,from,comm);
      DEBUG_OUT(LVL2,(stdout,"%s: leaving DDI_Recv.\n",DDI_Id()))
   }


/* -------------------------------------------------------------------- *\
   DDI_Comm_recv(bufer,size,from,comm)
   ===================================
   [OUT] buff - address in which incoming data is placed.
   [IN]  size - size of incoming data stream.
   [IN]  from - rank of process sending the data stream.
   [IN]  comm - DDI communicator Id.
   
   Used for point-to-point communications and large internal data
   tranfers to/from data servers.
\* -------------------------------------------------------------------- */
   void DDI_Recv_comm(void *buff,size_t size,int from,int comm_id) {
      const DDI_Comm *comm = (DDI_Comm *) Comm_find(comm_id);
      DEBUG_OUT(LVL1,(stdout,"%s: entering DDI_Comm_recv.\n",DDI_Id()))
      Comm_recv(buff,size,from,comm);
      DEBUG_OUT(LVL2,(stdout,"%s: leaving DDI_Comm_recv.\n",DDI_Id()))
   }


/* -------------------------------------------------------------------- *\
   Used to send a point-to-point message of a capped size.
\* -------------------------------------------------------------------- */
   void Comm_recv(void *buff,size_t size,int from,const DDI_Comm *comm) {

      char *buffer = (char *) buff;
      size_t working_size;
      size_t remaining = size;
      size_t max_recv = MAX_PTP_MSG_SIZE;

      DEBUG_OUT(LVL4,(stdout,"%s: entering Comm_recv(%x,%li,%i).\n",DDI_Id(),buff,size,from))

    # if defined DDI_CHECK_ARGS
      if(from < 0 || from > gv(ddi_base_comm).np*2) {
         fprintf(stdout,"%s: Comm_recv cannot recv data from process %i.\n",DDI_Id(),from);
         Fatal_error(911);
      }
    # endif
      
      while(remaining) {

         working_size = max_recv;
         if(remaining < working_size) working_size = remaining;

         Comm_recv_(buffer,working_size,from,comm);

         remaining -= working_size;
         buffer    += working_size;

      }

      DEBUG_OUT(LVL4,(stdout,"%s: leaving Comm_recv.\n",DDI_Id()))
   }


/* -------------------------------------------------------------------- *\
   Used to send a point-to-point message.
\* -------------------------------------------------------------------- */
   void Comm_recv_(void *buff,size_t size,int from,const DDI_Comm *comm) {

    # if defined DDI_MPI
      MPI_Status status;
      int size_i4 = (int) size;
    # if defined DDI_SOC && defined USE_TCP_LOCAL
      int global_cp;
      const DDI_Comm *comm_world = (DDI_Comm *) Comm_find(DDI_COMM_WORLD);
    # endif
    # endif

    # if defined SOC_SYNC
      char ack = 37;
    # endif

   /* ----------------------------------------- *\
    * Determine any mapping from communicators
   \* ----------------------------------------- */
      int global_pid = from;
      if(from < comm->np) global_pid = comm->global_pid[from];

    # if defined DDI_MPI
    # if defined DDI_SOC && defined USE_TCP_LOCAL
   /* -------------------------------------------- *\
    * Use TCP for local intra-node communication.
    * - skip MPI_Recv if incoming message is local
   \* -------------------------------------------- */
      global_cp = global_pid;
      if(global_cp >= comm_world->np) global_cp -= comm_world->np;
      if(comm_world->local_nid[global_cp] != comm_world->my) {
    # endif
         DEBUG_OUT(LVL8,(stdout,"%s: comm_recv_ - mpi_recv(buff,%li,MPI_BYTE,%i,0,%i,status)\n",
                         DDI_Id(),size,global_pid,comm->world_comm))
         if ( (size > DDI_LG_MSG_THRESHOLD) && (DDI_LG_MSG_FACTOR >1) ) {
            MPI_Datatype ddi_msg_type;
            MPI_Type_contiguous(DDI_LG_MSG_FACTOR, MPI_BYTE, &ddi_msg_type);
            MPI_Type_commit(&ddi_msg_type);
            size_t sizerecv;
            sizerecv = size/DDI_LG_MSG_FACTOR;
            //check message count isn't too big
            if (sizerecv > INT_MAX){
              fprintf(stdout,"WARNING!!!: Size of array to be sent sizerecv=%ld through MPI_Isend exceeded INT_MAX=%i\n",
                        sizerecv, INT_MAX);
              fprintf(stdout,"This situation is not supported by integer datatype for count of send buffer according to MPI standards,\n");
              fprintf(stdout,"therefore this job has been stopped.\n");
              fprintf(stdout,"You MUST either reduce the MEMORY= or MWORDS= value,\n");
              fprintf(stdout,"         or recompile the ddi source code with increased value of DDI_MSG_FACTOR.\n");
              fprintf(stdout,"This binary was compiled with DDI_MSG_FACTOR=%i\n",DDI_LG_MSG_FACTOR);
              Fatal_error(911);
            }
            if (size % DDI_LG_MSG_FACTOR != 0){
              // create new data type that is a grouping of first type
              //              //    for the first part of the message...this is the largest
              MPI_Datatype first_ddi_msg_type;
              MPI_Type_contiguous(sizerecv, ddi_msg_type, &first_ddi_msg_type);
              MPI_Type_commit(&first_ddi_msg_type);
              // create a new data type for the remainder of the message
              //              //     or the tail
              int remainder = size - sizerecv*DDI_LG_MSG_FACTOR ;
              MPI_Datatype remainder_type;
              MPI_Type_contiguous(remainder, MPI_BYTE, &remainder_type);
              MPI_Type_commit(&remainder_type);
              // create a new data type that is the first and remainder combined
              MPI_Datatype full_msg_type;
              int lengths[2] = { 1, 1 };
              const MPI_Aint displacements[2] = { 0, sizerecv*DDI_LG_MSG_FACTOR*sizeof(char)  };
              MPI_Datatype types[2] = { first_ddi_msg_type, remainder_type };
              MPI_Type_create_struct(2, lengths, displacements, types, &full_msg_type);
              MPI_Type_commit(&full_msg_type);
              MPI_Recv(buff,1,full_msg_type,global_pid,0,comm->world_comm,&status);
              // free all the new types created
              MPI_Type_free (&first_ddi_msg_type);
              MPI_Type_free (&remainder_type);
              MPI_Type_free (&full_msg_type);
            }else{
              MPI_Recv(buff,sizerecv,ddi_msg_type,global_pid,0,comm->world_comm,&status);
            }
            MPI_Type_free (&ddi_msg_type);
         }else{
            MPI_Recv(buff,size_i4,MPI_BYTE,global_pid,0,comm->world_comm,&status);
         }
         return;
    # if defined DDI_SOC && defined USE_TCP_LOCAL
      }
    # endif
    # endif

      
    # if defined DDI_SOC
      DEBUG_OUT(LVL8,(stdout,"%s: comm_recv_ - tcp_recv\n",DDI_Id()))
      Recv(gv(sockets)[global_pid],buff,size,0);
    # if defined SOC_SYNC
      Send(gv(sockets)[global_pid],&ack,1,0);
    # endif
    # endif

    }
      

/* -------------------------------------------------------------------- *\
\* -------------------------------------------------------------------- */
   void Comm_recv_list(void *buff,size_t size,int from,const DDI_List *list) {

      char *buffer = (char *) buff;
      size_t working_size;
      size_t remaining = size;

      DEBUG_OUT(LVL6,(stdout,"%s: entering Comm_recv_list.\n",DDI_Id()))

    # if defined DDI_CHECK_ARGS
      if(from < 0 || from > list->np) {
         fprintf(stdout,"%s: Comm_recv_list cannot recv data from process %i.\n",DDI_Id(),from);
         Fatal_error(911);
      }
    # endif
      
      while(remaining) {

         working_size = MAX_PTP_MSG_SIZE;
         if(remaining < working_size) working_size = remaining;

         Comm_recv_list_(buffer,working_size,from,list);

         remaining -= working_size;
         buffer    += working_size;

      }

      DEBUG_OUT(LVL6,(stdout,"%s: leaving Comm_recv_list.\n",DDI_Id()))
   }


/* -------------------------------------------------------------------- *\
\* -------------------------------------------------------------------- */
   void Comm_recv_list_(void *buff,size_t size,int from,const DDI_List *list) {

    # if defined DDI_SOC
    # if defined SOC_SYNC
      char ack = 31;
    # endif
      int global_pid = list->pids[from];
      Recv(gv(sockets)[global_pid],buff,size,0);
    # if defined SOC_SYNC
      Send(gv(sockets)[global_pid],&ack,1,0);
    # endif
    # endif

    # if defined DDI_MPI
      MPI_Status status;

      if ( (size > DDI_LG_MSG_THRESHOLD) && (DDI_LG_MSG_FACTOR >1) ) {
         MPI_Datatype ddi_msg_type;
         MPI_Type_contiguous(DDI_LG_MSG_FACTOR, MPI_BYTE, &ddi_msg_type);
         MPI_Type_commit(&ddi_msg_type);
         size_t sizerecv;
         sizerecv = size/DDI_LG_MSG_FACTOR;
         //check message count isn't too big
         if (sizerecv > INT_MAX){
           fprintf(stdout,"WARNING!!!: Size of array to be sent sizerecv=%ld through MPI_Isend exceeded INT_MAX=%i\n",
                        sizerecv, INT_MAX);
           fprintf(stdout,"This situation is not supported by integer datatype for count of send buffer according to MPI standards,\n");
           fprintf(stdout,"therefore this job has been stopped.\n");
           fprintf(stdout,"You MUST either reduce the MEMORY= or MWORDS= value,\n");
           fprintf(stdout,"         or recompile the ddi source code with increased value of DDI_MSG_FACTOR.\n");
           fprintf(stdout,"This binary was compiled with DDI_MSG_FACTOR=%i\n",DDI_LG_MSG_FACTOR);
           Fatal_error(911);
         }
         if (size % DDI_LG_MSG_FACTOR != 0){
             // create new data type that is a grouping of first type
             //              //    for the first part of the message...this is the largest
             MPI_Datatype first_ddi_msg_type;
             MPI_Type_contiguous(sizerecv, ddi_msg_type, &first_ddi_msg_type);
             MPI_Type_commit(&first_ddi_msg_type);
             // create a new data type for the remainder of the message
             //              //     or the tail
             int remainder = size - sizerecv*DDI_LG_MSG_FACTOR ;
             MPI_Datatype remainder_type;
             MPI_Type_contiguous(remainder, MPI_BYTE, &remainder_type);
             MPI_Type_commit(&remainder_type);
             // create a new data type that is the first and remainder combined
             MPI_Datatype full_msg_type;
             int lengths[2] = { 1, 1 };
             const MPI_Aint displacements[2] = { 0, sizerecv*DDI_LG_MSG_FACTOR*sizeof(char)  };
             MPI_Datatype types[2] = { first_ddi_msg_type, remainder_type };
             MPI_Type_create_struct(2, lengths, displacements, types, &full_msg_type);
             MPI_Type_commit(&full_msg_type);
             MPI_Recv(buff,1,full_msg_type,from,0,list->comm,&status);
             // free all the new types created
             MPI_Type_free (&first_ddi_msg_type);
             MPI_Type_free (&remainder_type);
             MPI_Type_free (&full_msg_type);
         }else{
             MPI_Recv(buff,sizerecv,ddi_msg_type,from,0,list->comm,&status);
         }
         MPI_Type_free (&ddi_msg_type);
     }else{
         MPI_Recv(buff,size,MPI_BYTE,from,0,list->comm,&status);

     }

    # endif

   }

