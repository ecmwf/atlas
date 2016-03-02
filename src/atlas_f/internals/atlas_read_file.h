
#define Char char

extern "C"
{
int atlas__read_file (const char* path, Char* &content, int &size);
}
#undef Char
