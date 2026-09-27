/* -------------------------------------------------------------------- *\
 * Distributed Data Interface
 * ==========================
 * 
 * Global objects
 * 
 * Author: Ryan M. Olson
 * CVS $Id: ddi_recv.c,v 1.1.1.1 2007/05/26 01:42:29 andrey Exp $
\* -------------------------------------------------------------------- */
 # include "ddi_base.h"
 

/* ------------------------------------------------------------------------- *\
   DDI_Recv_request(buff,from)
   ===========================
   [OUT] buff - address in which incoming data is passed.
   [IN]  from - rank of process sending data stream.
   
   Used to receive a data request sent by DDI_Send_request.  This subroutine
   is only called from a data server.
\* ------------------------------------------------------------------------- */
   void DDI_Recv_request(void *buff,int *from) {

    # if defined DDI_SOC
      char ack=37;
      size_t size = sizeof(DDI_Patch);
      DDI_Recvany(buff,size,from);
      Send(gv(sockets)[*from],&ack,1,0);
    # endif

   /* ---------------------------------- *\
      The stand-alone MPI version of DDI
   \* ---------------------------------- */

   /*
      This routine is the one that is responsible for the poor
      performance of a 100% MPI-1 model.  The "receive from anywhere"
      option caused by MPI_ANY_SOURCE typically is implemented by
      repeated checking (polling) on all open MPI processes.  This
      can sometimes be influenced by looking for options to control
      the polling mechanism, at the "mpirun" time.  However, a more
      practical solution is to use the "mixed" model, which uses
      our TCP/IP code for handling the small control messages, which
      do no polling at all.  If your adapter allows TCP/IP to coexist
      with the MPI-1, by all means try "mixed" over "mpi" for DDI.

      Here's a short note from Ryan, 
         Having talking to some MPI people at past SC conferences,
         it seems that the Irecv/wait combination has the best 
         chance to be implemented in a blocking fashion.  However, 
         I never worked out all the details to make it happen.

      What I think this means is that you set up an IRECV on every
      possible source, and then loop around "polling" with the wait
      yourself, perhaps with a delay from "sleep" added before the
      loop repeats itself.

    */

    # if defined DDI_MPI && !defined DDI_SOC & !defined DDI_LAPI
      const DDI_Comm *comm = (const DDI_Comm *) Comm_find(DDI_WORKING_COMM);
      MPI_Status status;
      size_t size = sizeof(DDI_Patch);
      DEBUG_OUT(LVL4,(stdout,"%s: call mpi_recv from any source.\n",DDI_Id()))


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
             //    for the first part of the message...this is the largest
             MPI_Datatype first_ddi_msg_type;
             MPI_Type_contiguous(sizerecv, ddi_msg_type, &first_ddi_msg_type);
             MPI_Type_commit(&first_ddi_msg_type);
             // create a new data type for the remainder of the message
             //     or the tail
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
             MPI_Recv(buff,1,full_msg_type,MPI_ANY_SOURCE,0,comm->world_comm,&status);
             // free all the new types created
             MPI_Type_free (&first_ddi_msg_type);
             MPI_Type_free (&remainder_type);
             MPI_Type_free (&full_msg_type);
         }else{
              MPI_Recv(buff,sizerecv,ddi_msg_type,MPI_ANY_SOURCE,0,comm->world_comm,&status);
         }
         MPI_Type_free (&ddi_msg_type);
     }else{
          MPI_Recv(buff,size,MPI_BYTE,MPI_ANY_SOURCE,0,comm->world_comm,&status);
     }
      *from = status.MPI_SOURCE;
      DEBUG_OUT(LVL4,(stdout,"%s: received request from %i.\n",DDI_Id(),*from))
    # endif

   }
