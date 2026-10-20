/* -------------------------------------------------------------------- *\
 * Distributed Data Interface
 * ==========================
 * 
 * Subroutines associated with the send operation using a communicator.
 *
 * Author: Ryan M. Olson
 * 10 Jun 09 - RMO - allow for one data server/node on Cray
 * 18 Aug 10 - MWS - correct upper rank calculation in 1st arg checking
\* -------------------------------------------------------------------- */
 # include "ddi_base.h"


/* -------------------------------------------------------------------- *\
   DDI_Send(bufer,size,from)
   =========================
   [IN] buff - address in which incoming data is placed.
   [IN] size - size of incoming data stream.
   [IN] to   - rank of process sending the data stream.
   
   Used for point-to-point communications and large internal data
   tranfers to/from data servers.  Use the working communicator.
\* -------------------------------------------------------------------- */
   void DDI_Send(void* buff,size_t size,int to) {
      const DDI_Comm *comm = (DDI_Comm *) Comm_find(DDI_WORKING_COMM);
      DEBUG_OUT(LVL1,(stdout,"%s: entering DDI_Send.\n",DDI_Id()))
      Comm_send(buff,size,to,comm);
      DEBUG_OUT(LVL2,(stdout,"%s: leaving DDI_Send.\n",DDI_Id()))
   }


/* -------------------------------------------------------------------- *\
   DDI_Comm_send(bufer,size,from,comm)
   ===================================
   [IN] buff - address in which incoming data is placed.
   [IN] size - size of incoming data stream.
   [IN] to   - rank of process sending the data stream.
   [IN] comm - DDI communicator Id.
   
   Used for point-to-point communications and large internal data
   tranfers to/from data servers.
\* -------------------------------------------------------------------- */
   void DDI_Send_comm(void* buff,size_t size,int to,int comm_id) {
      const DDI_Comm *comm = (DDI_Comm *) Comm_find(comm_id);
      DEBUG_OUT(LVL1,(stdout,"%s: entering DDI_Comm_send.\n",DDI_Id()))
      Comm_send(buff,size,to,comm);
      DEBUG_OUT(LVL2,(stdout,"%s: leaving DDI_Comm_send.\n",DDI_Id()))
   }


/* -------------------------------------------------------------------------------------- *\
   Used to send DDI messages of a capped size.
\* -------------------------------------------------------------------------------------- */
   void Comm_send(void *buff,size_t size,int to,const DDI_Comm *comm) {

      char ack = 37;
      char *buffer = (char *) buff;
      size_t working_size;
      size_t remaining = size;
      size_t max_send = MAX_PTP_MSG_SIZE;

      DEBUG_OUT(LVL4,(stdout,"%s: entering Comm_send(%x,%li,%i).\n",DDI_Id(),buff,size,to))

    # if defined DDI_CHECK_ARGS
      if(to < 0 || to > gv(ddi_base_comm).np*2) {
         fprintf(stdout,"%s: Comm_send cannot send data to process %i.\n",
                            DDI_Id(),to);
         Fatal_error(911);
      }
    # endif
           
      while(remaining) {

         working_size = max_send;
         if(remaining < working_size) working_size = remaining;

         Comm_send_(buffer,working_size,to,comm);

         remaining -= working_size;
         buffer    += working_size;

      }

      DEBUG_OUT(LVL4,(stdout,"%s: leaving Comm_send.\n",DDI_Id()))
   }


/* -------------------------------------------------------------------------------------- *\
   Used to send DDI messages.
\* -------------------------------------------------------------------------------------- */
   void Comm_send_(void* buff,size_t size,int to,const DDI_Comm *comm) {

    # if defined DDI_MPI
      int size_i4 = (int) size;
    # if defined DDI_SOC && defined USE_TCP_LOCAL
      int global_cp;
      const DDI_Comm *comm_world = (DDI_Comm *) Comm_find(DDI_COMM_WORLD);
    # endif
    # endif

    # if defined SOC_SYNC
      char ack = 3;
    # endif

   /* ----------------------------------------- *\
    * Determine any mapping from communicators
   \* ----------------------------------------- */
      int global_pid = to;
      if(to < comm->np) global_pid = comm->global_pid[to];

    # if defined DDI_MPI
    # if defined DDI_SOC && defined USE_TCP_LOCAL
   /* -------------------------------------------- *\
    * Use TCP for local intra-node communication.
    * - skip MPI_Send if outgoing message is local
   \* -------------------------------------------- */
      global_cp = global_pid;
      if(global_cp >= comm_world->np) global_cp -= comm_world->np;
      if(comm_world->local_nid[global_cp] != comm_world->my) {
    # endif
         DEBUG_OUT(LVL8,(stdout,
            "%s: comm_send_ - mpi_ssend(buff,%li,MPI_BYTE,%i,0,%i)\n",
            DDI_Id(),size,global_pid,comm->world_comm))

         if ( (size > DDI_LG_MSG_THRESHOLD) && (DDI_LG_MSG_FACTOR >1) ) {
            // Create the base datatype....this can probably be moved to ddi init
            MPI_Datatype ddi_msg_type;
            MPI_Type_contiguous(DDI_LG_MSG_FACTOR, MPI_BYTE, &ddi_msg_type);
            MPI_Type_commit(&ddi_msg_type);
            size_t sizesend;
            sizesend = size/DDI_LG_MSG_FACTOR;
            //check message count isn't too big
            if (sizesend > INT_MAX){
                fprintf(stdout,"WARNING!!!: Size of array to be sent sizesend=%ld through MPI_Isend exceeded INT_MAX=%d\n",sizesend, INT_MAX);
                fprintf(stdout,"This situation is not supported by integer datatype for count of send buffer according to MPI standards,\n");
                fprintf(stdout,"therefore this job has been stopped.\n");
                fprintf(stdout,"You MUST either reduce the MEMORY= or MWORDS= value,\n");
                fprintf(stdout,"         or recompile the ddi source code with increased value of DDI_MSG_FACTOR.\n");
                fprintf(stdout,"This binary was compiled with DDI_MSG_FACTOR=%i\n",DDI_LG_MSG_FACTOR);
                Fatal_error(911);
            }
            if (size % DDI_LG_MSG_FACTOR != 0){
                // create new data type that is a grouping of first type
                //    for the first part of the message...this is the largest
                MPI_Datatype first_ddi_msg_type;
                MPI_Type_contiguous(sizesend, ddi_msg_type, &first_ddi_msg_type);
                MPI_Type_commit(&first_ddi_msg_type);
                // create a new data type for the remainder of the message
                //     or the tail
                int remainder = size - sizesend*DDI_LG_MSG_FACTOR ;
                MPI_Datatype remainder_type;
                MPI_Type_contiguous(remainder, MPI_BYTE, &remainder_type);
                MPI_Type_commit(&remainder_type);
                // create a new data type that is the first and remainder combined
                MPI_Datatype full_msg_type;
                int lengths[2] = { 1, 1 };
                const MPI_Aint displacements[2] = { 0, sizesend*DDI_LG_MSG_FACTOR*sizeof(char)  };
                MPI_Datatype types[2] = { first_ddi_msg_type, remainder_type };
                MPI_Type_create_struct(2, lengths, displacements, types, &full_msg_type);
                MPI_Type_commit(&full_msg_type);
                MPI_Send(buff,1,full_msg_type,global_pid,0,comm->world_comm);
                // free all the new types created
                MPI_Type_free (&first_ddi_msg_type);
                MPI_Type_free (&remainder_type);
                MPI_Type_free (&full_msg_type);
            }else{
                MPI_Send(buff,sizesend,ddi_msg_type,global_pid,0,comm->world_comm);
            }
            MPI_Type_free (&ddi_msg_type);
         }else{
            MPI_Send(buff,size_i4,MPI_BYTE,global_pid,0,comm->world_comm);
         }

         return;
    # if defined DDI_SOC && defined USE_TCP_LOCAL
      }
    # endif
    # endif


    # if defined DDI_SOC
      DEBUG_OUT(LVL8,(stdout,"%s: comm_send_ - tcp_send\n",DDI_Id()))
      Send(gv(sockets)[global_pid],buff,size,0);
    # if defined SOC_SYNC
      Recv(gv(sockets)[global_pid],&ack,1,0);
      if(ack != 37) {
         fprintf(stdout,"%s: SOC_ERROR (ack=%i).\n",DDI_Id(),(int) ack);
         Fatal_error(911);
      }
    # endif
    # endif

   }


/* -------------------------------------------------------------------------------------- *\
\* -------------------------------------------------------------------------------------- */
   void Comm_send_list(void *buff,size_t size,int to,const DDI_List *list) {

      char *buffer = (char *) buff;
      size_t working_size;
      size_t remaining = size;

      DEBUG_OUT(LVL6,(stdout,"%s: entering Comm_send_list.\n",DDI_Id()))

    # if defined DDI_CHECK_ARGS
      if(to < 0 || to > list->np) {
         fprintf(stdout,"%s: Comm_send_list cannot send data to process %i.\n",DDI_Id(),to);
         Fatal_error(911);
      }
    # endif
           
      while(remaining) {

         working_size = MAX_PTP_MSG_SIZE;
         if(remaining < working_size) working_size = remaining;

         Comm_send_list_(buffer,working_size,to,list);

         remaining -= working_size;
         buffer    += working_size;

      }

      DEBUG_OUT(LVL6,(stdout,"%s: leaving Comm_send_list.\n",DDI_Id()))
   }


/* -------------------------------------------------------------------------------------- *\
\* -------------------------------------------------------------------------------------- */
   void Comm_send_list_(void* buff,size_t size,int to,const DDI_List *list) {

    # if defined DDI_SOC
    # if defined SOC_SYNC
      char ack = 3;
    # endif
      int global_pid = list->pids[to];
      Send(gv(sockets)[global_pid],buff,size,0);
    # if defined SOC_SYNC
      Recv(gv(sockets)[global_pid],&ack,1,0);
      if(ack != 31) {
         fprintf(stdout,"%s: SOC_ERROR (ack=%i).\n",DDI_Id(),(int) ack);
         Fatal_error(911);
      }
    # endif
    # endif

    # if defined DDI_MPI
      if ( (size > DDI_LG_MSG_THRESHOLD) && (DDI_LG_MSG_FACTOR >1) ) {
         // Create the base datatype....this can probably be moved to ddi init
         MPI_Datatype ddi_msg_type;
         MPI_Type_contiguous(DDI_LG_MSG_FACTOR, MPI_BYTE, &ddi_msg_type);
         MPI_Type_commit(&ddi_msg_type);
         size_t sizesend;
         sizesend = size/DDI_LG_MSG_FACTOR;
         //check message count isn't too big
         if (sizesend > INT_MAX){
             fprintf(stdout,"WARNING!!!: Size of array to be sent sizesend=%ld through MPI_Isend exceeded INT_MAX=%d\n",sizesend, INT_MAX);
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
             MPI_Type_contiguous(sizesend, ddi_msg_type, &first_ddi_msg_type);
             MPI_Type_commit(&first_ddi_msg_type);
             // create a new data type for the remainder of the message
             //              //     or the tail
             int remainder = size - sizesend*DDI_LG_MSG_FACTOR ;
             MPI_Datatype remainder_type;
             MPI_Type_contiguous(remainder, MPI_BYTE, &remainder_type);
             MPI_Type_commit(&remainder_type);
             // create a new data type that is the first and remainder combined
             MPI_Datatype full_msg_type;
             int lengths[2] = { 1, 1 };
             const MPI_Aint displacements[2] = { 0, sizesend*DDI_LG_MSG_FACTOR*sizeof(char)  };
             MPI_Datatype types[2] = { first_ddi_msg_type, remainder_type };
             MPI_Type_create_struct(2, lengths, displacements, types, &full_msg_type);
             MPI_Type_commit(&full_msg_type);
             MPI_Ssend(buff,1,full_msg_type,to,0,list->comm);
             // free all then new types created
             MPI_Type_free (&first_ddi_msg_type);
             MPI_Type_free (&remainder_type);
             MPI_Type_free (&full_msg_type);
         }else{
             MPI_Ssend(buff,sizesend,ddi_msg_type,to,0,list->comm);
         }
         MPI_Type_free (&ddi_msg_type);
      }else{
         MPI_Ssend(buff,size,MPI_BYTE,to,0,list->comm);
      }
    # endif

   }

