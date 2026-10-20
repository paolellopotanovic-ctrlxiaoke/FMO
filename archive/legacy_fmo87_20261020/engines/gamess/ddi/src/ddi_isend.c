/* ------------------------------------------------------------------ *\
   Subroutine DDI_ISend(BUFF,SIZE,TYPE,TO,TAG,REQ)
   ===============================================
   [IN]  BUFF - Buffer to be sent.
   [IN]  SIZE - Size of buffer in bytes.
   [IN]  TO   - Rank of processor (within current scope) to whom the
                message is to be sent.
   
   Author: Ryan M. Olson
   CVS $Id: ddi_isend.c,v 1.1.1.1 2007/05/26 01:42:30 andrey Exp $
\* ------------------------------------------------------------------ */
 # include "ddi_base.h"

   static void *DDI_ISend_thread(void *);
 
   void DDI_ISend(void *buffer,size_t size,int to,int *req_val) {

      int global_pid;
      DDI_Request *req = &gv(isend_req);
      const DDI_Comm *comm = (const DDI_Comm *) Comm_find(DDI_WORKING_COMM);
    
    # if defined DDI_CHECK_ARGS
      if(to < 0 || to > comm->np) {
         fprintf(stdout,"%s: can not send to ddi process %i.\n",DDI_Id(),to);
         Fatal_error(911);
      }
    # endif

      
    # if defined DDI_SOC && !defined DDI_MPI
      pthread_attr_t thread_attr;
      pthread_attr_init(&thread_attr);
      pthread_attr_setscope(&thread_attr,PTHREAD_SCOPE_SYSTEM);
      req->to     = to;
      req->size   = size;
      req->buffer = buffer;
      if(pthread_create(&req->hnd_thread,&thread_attr,DDI_ISend_thread,req) == -1) {
         fprintf(stderr,"%s: pthread_create failed in DDI_ISend.\n",DDI_Id());
         Fatal_error(911);
      }
    # endif
     
    # if defined DDI_MPI
      //check if size is larger that the ddi large message threshold
      //  and that the DDI_LG_MSG_FACTOR has been set, default is 1
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

         // check if size was evenly divisible by DDI_LG_MSG_FACTOR and if not create structure for remainder
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

             MPI_Isend(buffer,1,full_msg_type,to,1,comm->compute_comm,req);
             // free all then new types created
             MPI_Type_free (&first_ddi_msg_type);
             MPI_Type_free (&remainder_type);
             MPI_Type_free (&full_msg_type);

         }else{
             MPI_Isend(buffer,sizesend,ddi_msg_type,to,1,comm->compute_comm,req);
         }
         MPI_Type_free (&ddi_msg_type);

      }else{
         MPI_Isend(buffer,size,MPI_BYTE,to,1,comm->compute_comm,req);
      }
    # endif
     
      ULTRA_DEBUG((stdout,"%s: non-blocking send to %i issued.\n",DDI_Id(),to))

      *req_val = 1;
   }


 # ifndef DDI_MPI  
   static void *DDI_ISend_thread(void *myarg) {
      DDI_Request *req = (DDI_Request *) myarg;
      DDI_Send(req->buffer,req->size,req->to);
      ULTRA_DEBUG((stdout,"%s: isend_thread finished.\n",DDI_Id()))
      return NULL;
   }  
 # endif
 
