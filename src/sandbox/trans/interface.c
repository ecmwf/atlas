#include <stdio.h> /* printf */
#include <stdlib.h> /* malloc, free */

#include "transi/trans.h"

/*---------------------------------------------------------------------------------*/

/* in  "transi/trans.h" */

/// @returns NULL if not in cache
typedef void* (*open_proc)(void* ctxt, const char* key, int mode);

typedef size_t (*write_proc)(void* ctxt, const void* buffer, size_t size, void* data);
typedef size_t (*read_proc) (void* ctxt, void* buffer, size_t size, void* data);
typedef int (*close_proc)(void* ctxt, void* data);

typedef struct cache_handlers
{
	void* context;
	open_proc open;
	read_proc read;
	write_proc write;
	close_proc close;
} cache_handlers;

static cache_handlers* global_handlers;

// null means no caching
void trans_register_cache_handlers(cache_handlers* handlers)
{
	global_handlers = handlers;
}

void* trans_cached_open_(const char* key, int* mode/*, fortint* err*/)
{

}
size_t trans_cached_read_(const char* key, int* mode/*, fortint* err*/)
{

}

void compute_()
{
	printf("Computing\n");
}

void trans_compute_or_cache(struct Trans_t* trans)
{
	if(global_handlers)
	{
		void* ctxt = global_handlers->context;

		/* TODO: based on trans handle create the key */
		char key[256] = "T1279_rgg.N640";

		/* TODO: based on trans handle compute & allocate appropriate size for buffer */
		char buffer[1024];
		size_t sz = 1024;

		int mode = 0; // READ

		void* cache_entry = global_handlers->open(ctxt, key, mode);

		if( cache_entry )
		{
			global_handlers->read(ctxt, buffer, sz, cache_entry);
		}
		else
		{
			compute_();

			/* TODO: copy / fill buffers from trans handle */

			mode = 1; // WRITE
			global_handlers->open(ctxt, key, mode);
			global_handlers->write(ctxt, buffer, sz, cache_entry);
		}


		// ASSERT( cache_entry )

		global_handlers->close(ctxt, cache_entry);
	}
	else
		compute_();
}

/*---------------------------------------------------------------------------------*/

void* my_open(void* ctxt, const char* key, int mode)
{
	return fopen(key, mode ? "w" : "r");
}

size_t my_write(void* ctxt, const void* buffer, size_t size, void* data)
{
	FILE* f = (FILE*) data;
	return fwrite(buffer,1,size,f);
}

size_t my_read (void* ctxt, void* buffer, size_t size, void* data)
{
	FILE* f = (FILE*) data;
	return fread(buffer,1,size,f);
}

int my_close(void* ctxt, void* data)
{
	FILE* f = (FILE*) data;
	return fclose(f);
}

struct cache_handlers my_handlers = { NULL, &my_open, &my_read, &my_write, &my_close };

/*---------------------------------------------------------------------------------*/

void read_grid(struct Trans_t* trans)
{
  int i;
  int T159[] = {
     18,  25,  36,  40,  45,  50,  60,  64,  72,  72,
     80,  90,  96, 100, 108, 120, 120, 125, 135, 144,
    144, 150, 160, 160, 180, 180, 180, 192, 192, 200,
    200, 216, 216, 216, 225, 225, 240, 240, 240, 243,
    250, 250, 256, 270, 270, 270, 288, 288, 288, 288,
    288, 288, 300, 300, 300, 300, 320, 320, 320, 320,
    320, 320, 320, 320, 320, 320, 320, 320, 320, 320,
    320, 320, 320, 320, 320, 320, 320, 320, 320, 320,
    320, 320, 320, 320, 320, 320, 320, 320, 320, 320,
    320, 320, 320, 320, 320, 320, 320, 320, 320, 320,
    320, 320, 320, 320, 300, 300, 300, 300, 288, 288,
    288, 288, 288, 288, 270, 270, 270, 256, 250, 250,
    243, 240, 240, 240, 225, 225, 216, 216, 216, 200,
    200, 192, 192, 180, 180, 180, 160, 160, 150, 144,
    144, 135, 125, 120, 120, 108, 100,  96,  90,  80,
     72,  72,  64,  60,  50,  45,  40,  36,  25,  18,
  };
  trans->ndgl  = sizeof(T159)/sizeof(int);
  trans->nloen = malloc( sizeof(T159) );
  for( i=0; i<trans->ndgl; i++)  trans->nloen[i] = T159[i];

  // Assume Linear Grid
  trans->nsmax=(2*trans->ndgl-1)/2;
}

int main ( int arc, char **argv )
{
  	trans_init();

  	trans_register_cache_handlers(&my_handlers);

  	struct Trans_t trans;
    trans_new(&trans);
  	read_grid(&trans);

  	// setup my parameters in trans handle
  	//

  	trans_setup(&trans); // this uses my_reader

	  	/* called inside trans(i)_setup */
  		int mode = 0;
	  	void* cache_data = trans_cached_open_("FOO",&mode); /* 0 - read, 1 - write */





  	trans_finalize();
}
