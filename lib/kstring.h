
// Taken from:
// https://github.com/lh3/bwa/blob/master/kstring.h

#ifndef KSTRING_H
#define KSTRING_H

#include <stdlib.h>
#include <string.h>
#include <stdarg.h>
#include <stdio.h>

#ifndef kroundup32
#define kroundup32(x) (--(x), (x)|=(x)>>1, (x)|=(x)>>2, (x)|=(x)>>4, (x)|=(x)>>8, (x)|=(x)>>16, ++(x))
#endif

#ifndef KSTRING_T
#define KSTRING_T kstring_t
typedef struct __kstring_t {
	size_t l, m;
	char *s;
} kstring_t;
#endif

static inline char *ks_str(kstring_t *s)
{
	return s->s;
}

static inline size_t ks_len(kstring_t *s)
{
	return s->l;
}

static inline void ks_resize(kstring_t *s, size_t size)
{
	if (s->m < size) {
		s->m = size;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
}

static inline int kputsn(const char *p, int l, kstring_t *s)
{
	if (s->l + l + 1 >= s->m) {
		s->m = s->l + l + 2;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
	memcpy(s->s + s->l, p, l);
	s->l += l;
	s->s[s->l] = 0;
	return l;
}

static inline int kputs(const char *p, kstring_t *s)
{
	return kputsn(p, strlen(p), s);
}

static inline int kputc(int c, kstring_t *s)
{
	if (s->l + 1 >= s->m) {
		s->m = s->l + 2;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
	s->s[s->l++] = c;
	s->s[s->l] = 0;
	return c;
}

static inline int kputw(int c, kstring_t *s)
{
	char buf[16];
	int l, x;
	if (c == 0) return kputc('0', s);
	for (l = 0, x = c < 0? -c : c; x > 0; x /= 10) buf[l++] = x%10 + '0';
	if (c < 0) buf[l++] = '-';
	if (s->l + l + 1 >= s->m) {
		s->m = s->l + l + 2;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
	for (x = l - 1; x >= 0; --x) s->s[s->l++] = buf[x];
	s->s[s->l] = 0;
	return 0;
}

static inline int kputuw(unsigned c, kstring_t *s)
{
	char buf[16];
	int l, i;
	unsigned x;
	if (c == 0) return kputc('0', s);
	for (l = 0, x = c; x > 0; x /= 10) buf[l++] = x%10 + '0';
	if (s->l + l + 1 >= s->m) {
		s->m = s->l + l + 2;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
	for (i = l - 1; i >= 0; --i) s->s[s->l++] = buf[i];
	s->s[s->l] = 0;
	return 0;
}

static inline int kputl(long c, kstring_t *s)
{
	char buf[32];
	long l, x;
	if (c == 0) return kputc('0', s);
	for (l = 0, x = c < 0? -c : c; x > 0; x /= 10) buf[l++] = x%10 + '0';
	if (c < 0) buf[l++] = '-';
	if (s->l + l + 1 >= s->m) {
		s->m = s->l + l + 2;
		kroundup32(s->m);
		s->s = (char*)realloc(s->s, s->m);
	}
	for (x = l - 1; x >= 0; --x) s->s[s->l++] = buf[x];
	s->s[s->l] = 0;
	return 0;
}

// Free memory associated with a kstring_t*
static inline void destroy_kstring(kstring_t* k) {
	if (k == NULL)
		return;
	if (ks_str(k) != NULL)
		free(ks_str(k)); 
	free(k);
}

// Added on 3 May 2024.
// Overwrites contents of s with l characters from p.
static inline int ks_overwriten(const char* p, int l, kstring_t *s) {
	destroy_kstring(s);
	s = calloc(1, sizeof(kstring_t));
	return kputsn(p, l, s);
}

// Added on 3 May 2024.
// Overwrites contents of s with p.
static inline void ks_overwrite(const char* p, kstring_t *s) {
	destroy_kstring(s);
	s = calloc(1, sizeof(kstring_t));
	kputs(p, s);
}

#endif