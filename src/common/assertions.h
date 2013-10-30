#ifndef common_assertions_h
#define common_assertions_h

#ifndef NDEBUG
#define assert(ASSERTION) \
if( .not. (ASSERTION) ) then; \
write(0,*) "Assertion failed: ASSERTION"; \
call abort; \
endif;
#else
#define assert(ASSERTION)
#endif

#ifndef NDEBUG
#define assert_msg(ASSERTION,msg) \
if( .not. (ASSERTION) ) then; \
write(0,*) "Assertion failed: "//trim(msg); \
call abort; \
endif;
#else
#define assert_msg(ASSERTION,msg)
#endif

#endif 


#define assert_always(ASSERTION) \
if( .not. (ASSERTION) ) then; \
write(0,*) "Assertion failed: ASSERTION"; \
call abort; \
endif;
