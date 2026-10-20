/* -------------------------------------------------------------------- *\
 * Distributed Data Interface
 * ==========================
 * 
 * Subroutines associated with the send operation.
 *
 * Author: Ryan M. Olson
 * 10 Jun 09 - RMO - update ddi_send_request args
\* -------------------------------------------------------------------- */
 # include "ddi_base.h"


/* --------------------------------------------------------- *\
   DDI_Send_request(buff,to,request)
   =================================
   [IN] buff    - address containing a DDI_Patch object.
   [IN] to      - target DDI process to send the data request.
   [IN] request - a structure of type DDI_Request
   
   This subroutine is only called by a compute process, and
   is used only to transmit a data request to a data server.
   It uses TCP/IP, MPI-1, or LAPI to issue the data request.
\* --------------------------------------------------------- */
   void DDI_Send_request(void *buff,int *to,DDI_Request *req) {
      DDI_Send_request_comm(buff,to,req,DDI_WORKING_COMM);
   }
   
   void DDI_Send_request_comm(void *buff,int *to,DDI_Request *req,int commid) {
      char ack;
      size_t size = sizeof(DDI_Patch);
      int i,np,me,nn,my;

      const DDI_Comm *comm = (const DDI_Comm *) Comm_find(commid);
    
      np = comm->np;
      me = comm->me;
      nn = comm->nn;
      my = comm->my; 

      *to = comm->global_dsid[*to];
/*      fprintf(stdout,"%s: sending request to global process %i.\n",DDI_Id(),*to);*/
      DEBUG_OUT(LVL3,(stdout,"%s: sending request to global process %i.\n",DDI_Id(),*to))

   /* ------------------------------------------------------------ *\
      Using TCP/IP sockets, this is always a synchronous operation
   \* ------------------------------------------------------------ */
    # if defined DDI_SOC
      Send(gv(sockets)[*to],buff,size,0);
      Recv(gv(sockets)[*to],&ack,1,0);
    # endif


   /* -------------------------------------------------------------- *\
      Using LAPI, this sends an active message to the target process
      causing an interrupt signal to be issued.  The target process
      now acts like a data server and handles the data request. This
      call is non-blocking, because the once the active message is
      sent, the target process is in control of the data and the 
      originating compute process only needs to wait until all the
      target process have finished and tells 'this' originating
      process that it can continue.  These are slightly different
      for get, put and accumulates.
   \* -------------------------------------------------------------- */
    # if defined DDI_LAPI
      DDI_Patch *patch = (DDI_Patch *) buff;
      uint tgt = gv(lapi_map)[*to];
      void *hdr_hndlr = (void *) gv(am_hndlr)[tgt];
      void *udata = NULL;
      ulong udata_len = 0;
      lapi_cntr_t *org_cntr = (lapi_cntr_t *) patch->cp_lapi_cntr;
      lapi_cntr_t *tgt_cntr = NULL;
      lapi_cntr_t *cmpl_cntr = NULL;
   
      if(LAPI_Amsend(gv(lapi_hnd),tgt,hdr_hndlr,buff,size,udata,udata_len,
                     tgt_cntr,org_cntr,cmpl_cntr) != LAPI_SUCCESS) {
          fprintf(stdout,"%s: lapi_amsend error in ddi_send_request.\n",DDI_Id());
          Fatal_error(911);
      }
    # endif


   /* ---------------------------------- *\
      The stand-alone MPI version of DDI
   \* ---------------------------------- */
    # if defined CRAY_MPI
      if(req == NULL) {
         MPI_Send(buff,size,MPI_BYTE,*to,37,comm->world_comm);
      } else {
         MPI_Isend(buff,size,MPI_BYTE,*to,37,comm->world_comm,req);
      }
      return;
    # endif
    # if defined DDI_MPI && !defined DDI_SOC && !defined DDI_LAPI
      DEBUG_OUT(LVL3,(stdout,"%s: calling mpi_ssend.\n",DDI_Id()))


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

             MPI_Ssend(buff,1,full_msg_type,*to,0,comm->world_comm);
             // free all then new types created
             MPI_Type_free (&first_ddi_msg_type);
             MPI_Type_free (&remainder_type);
             MPI_Type_free (&full_msg_type);
         }else{
             MPI_Ssend(buff,sizesend,ddi_msg_type,*to,0,comm->world_comm);
         }
         MPI_Type_free (&ddi_msg_type);
      }else{
         MPI_Ssend(buff,size,MPI_BYTE,*to,0,comm->world_comm);
      }

    # endif


   /* ------------------------------------------------------------- *\
      Reverse look up the group rank of the global rank to which to
      the data request is being sent.
   \* ------------------------------------------------------------- */
   /*
      if(gv(scope) == DDI_GROUP) {
         for(i=0; i<np; i++) if(DDI_Id_global_proc(i+np) == *to) *to = i+np;
      }
   */

      DEBUG_OUT(LVL3,(stdout,"%s: leaving ddi_send_request.\n",DDI_Id()))

   }

/* DDI_Send_request0 and DDI_Send_request_comm0 are clones of the two
   subroutines above; the only difference is the following line commented out:
   *to = comm->global_dsid[*to];
   It is needed when global_dsid is not available, and instead caller
   is smart enough to provide "to" to be the data server id.
*/

   void DDI_Send_request0(void *buff,int *to,DDI_Request *req) {
      DDI_Send_request_comm0(buff,to,req,DDI_WORKING_COMM);
   }
   
   void DDI_Send_request_comm0(void *buff,int *to,DDI_Request *req,int commid) {
      char ack;
      size_t size = sizeof(DDI_Patch);
      int i,np,me,nn,my;

      const DDI_Comm *comm = (const DDI_Comm *) Comm_find(commid);
    
      np = comm->np;
      me = comm->me;
      nn = comm->nn;
      my = comm->my; 

/*      *to = comm->global_dsid[*to]; */
/*      fprintf(stdout,"%s: sending request to global process %i.\n",DDI_Id(),*to);*/
      DEBUG_OUT(LVL3,(stdout,"%s: sending request to global process %i.\n",DDI_Id(),*to))

   /* ------------------------------------------------------------ *\
      Using TCP/IP sockets, this is always a synchronous operation
   \* ------------------------------------------------------------ */
    # if defined DDI_SOC
      Send(gv(sockets)[*to],buff,size,0);
      Recv(gv(sockets)[*to],&ack,1,0);
    # endif


   /* -------------------------------------------------------------- *\
      Using LAPI, this sends an active message to the target process
      causing an interrupt signal to be issued.  The target process
      now acts like a data server and handles the data request. This
      call is non-blocking, because the once the active message is
      sent, the target process is in control of the data and the 
      originating compute process only needs to wait until all the
      target process have finished and tells 'this' originating
      process that it can continue.  These are slightly different
      for get, put and accumulates.
   \* -------------------------------------------------------------- */
    # if defined DDI_LAPI
      DDI_Patch *patch = (DDI_Patch *) buff;
      uint tgt = gv(lapi_map)[*to];
      void *hdr_hndlr = (void *) gv(am_hndlr)[tgt];
      void *udata = NULL;
      ulong udata_len = 0;
      lapi_cntr_t *org_cntr = (lapi_cntr_t *) patch->cp_lapi_cntr;
      lapi_cntr_t *tgt_cntr = NULL;
      lapi_cntr_t *cmpl_cntr = NULL;
   
      if(LAPI_Amsend(gv(lapi_hnd),tgt,hdr_hndlr,buff,size,udata,udata_len,
                     tgt_cntr,org_cntr,cmpl_cntr) != LAPI_SUCCESS) {
          fprintf(stdout,"%s: lapi_amsend error in ddi_send_request.\n",DDI_Id());
          Fatal_error(911);
      }
    # endif


   /* ---------------------------------- *\
      The stand-alone MPI version of DDI
   \* ---------------------------------- */
    # if defined CRAY_MPI
      if(req == NULL) {
         MPI_Send(buff,size,MPI_BYTE,*to,37,comm->world_comm);
      } else {
         MPI_Isend(buff,size,MPI_BYTE,*to,37,comm->world_comm,req);
      }
      return;
    # endif
    # if defined DDI_MPI && !defined DDI_SOC && !defined DDI_LAPI
      DEBUG_OUT(LVL3,(stdout,"%s: calling mpi_ssend.\n",DDI_Id()))

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
             MPI_Ssend(buff,1,full_msg_type,*to,0,comm->world_comm);
             // free all then new types created
             MPI_Type_free (&first_ddi_msg_type);
             MPI_Type_free (&remainder_type);
             MPI_Type_free (&full_msg_type);
         }else{
             MPI_Ssend(buff,sizesend,ddi_msg_type,*to,0,comm->world_comm);
         }
         MPI_Type_free (&ddi_msg_type);
      }else{
         MPI_Ssend(buff,size,MPI_BYTE,*to,0,comm->world_comm);
      }

    # endif


   /* ------------------------------------------------------------- *\
      Reverse look up the group rank of the global rank to which to
      the data request is being sent.
   \* ------------------------------------------------------------- */
   /*
      if(gv(scope) == DDI_GROUP) {
         for(i=0; i<np; i++) if(DDI_Id_global_proc(i+np) == *to) *to = i+np;
      }
   */

      DEBUG_OUT(LVL3,(stdout,"%s: leaving ddi_send_request.\n",DDI_Id()))

   }

