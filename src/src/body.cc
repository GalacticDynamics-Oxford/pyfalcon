// -*- C++ -*-
////////////////////////////////////////////////////////////////////////////////
///
/// \file    src/public/lib/body.cc
///
/// \author  Walter Dehnen
///
/// \date    2000-2012
///
////////////////////////////////////////////////////////////////////////////////
//
// Copyright (C) 2000-2012  Walter Dehnen
//
// This program is free software; you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation; either version 2 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License along
// with this program; if not, write to the Free Software Foundation, Inc., 675
// Mass Ave, Cambridge, MA 02139, USA.
//
////////////////////////////////////////////////////////////////////////////////
#include <body.h>                                  // falcON::bodies etc
#include <iostream>                                // C++ basic I/O 
#include <fstream>                                 // C++ file I/O
#include <sstream>                                 // C++ string I/O
#include <iomanip>                                 // C++ I/O formating
#include <cstring>                                 // C++ strings 
//#include <public/nemo++.h>                         // utilities for NEMO I/O
//#include <public/bodyfunc.h>
#include <utils/numerics.h>
#include <utils/heap.h>

namespace falcON {
  using namespace WDutils;
}

using namespace falcON;
typedef long unsigned lu; // will convert size_t to this type in printf

falcON_TRAITS(falcON::bodies::block,"bodies::block");
//
// struct falcON::bodies::block
//
void bodies::block::clone(block*that)
{
  if(that == this) return;
  DebugInfo(3,"bodies::block::clone(): cloning block with %d [%d] %s\n",
	    that->NBOD,that->NALL,that->TYPE.name());
  if(this->TYPE != that->TYPE)
    falcON_THROW("bodies::block::clone(): bodytype mismatch ('%s' vs '%s')\n",
		 this->TYPE.name(), that->TYPE.name());
  for(fieldbit f; f; ++f) {
    this->del_field(f);
    this->set_data_void(f,that->data_void(f));
    that->set_data_void(f,0);
  }
  this->NALL       = that->NALL;
  this->NBOD       = that->NBOD;
  this->FIRST      = that->FIRST;
  this->LOCALFIRST = that->LOCALFIRST;
}
//
void bodies::block::reset_flags() const
{
  if(0 != DATA[fieldbit::f]) {
    if(TYPE.is_sph()) 
      for(unsigned n=0; n!=NALL; ++n)
	datum<fieldbit::f>(n) = flags::sph;
    else
    if(TYPE.is_sink()) 
      for(unsigned n=0; n!=NALL; ++n)
	datum<fieldbit::f>(n) = flags::sink;
    else
      for(unsigned n=0; n!=NALL; ++n)
	datum<fieldbit::f>(n) = flags::empty;
  }
}
//
void bodies::block::flag_all_as_active() const falcON_THROWING
{
  if(0 != DATA[fieldbit::f])
    for(unsigned n=0; n!=NALL; ++n)
      datum<fieldbit::f>(n).add(flags::active);
  else 
    falcON_THROW("in bodies::flag_all_as_active(): flags not supported");
}
//
void bodies::block::reset_data(fieldset b) const falcON_THROWING
{
#define RESETDATA(BIT,NAME)				\
  if(DATA[BIT] && b.contain(BIT) && NBOD)		\
    for(unsigned n=0; n!=NBOD; ++n)			\
      field_traits<BIT>::set_zero(datum<BIT>(n));
 DEF_NAMED(RESETDATA)
#undef RESETDATA
}
//
void bodies::block::add_field (fieldbit f) falcON_THROWING
{
  if(TYPE.allows(f) && 0 == DATA[value(f)] ) {
    DebugInfo(4,"bodies::block::add_field(): "
	      "allocating data for %s bodies: %u %c (%s)\n",
	      TYPE.name(),NALL,letter(f),fullname(f));
    set_data_void(f, falcON_NEW(char,NALL*falcON::size(f)));
    if(f == fieldbit::f) reset_flags();
  }
}
//
void bodies::block::del_field (fieldbit f) falcON_THROWING
{
  if(DATA[value(f)]) {
    DebugInfo(4,"bodies::block::del_field(): "
	      "de-allocating data for %s bodies: %c (%s)\n",
	      TYPE.name(),letter(f),fullname(f));
    falcON_DEL_A(static_cast<char*>(DATA[value(f)]));
  }
  set_data_void(f,0);
}
//
#if 0
void bodies::block::swap_bytes(fieldbit f) falcON_THROWING
{
  if(DATA[value(f)]) {
    DebugInfo(4,"bodies::block::swap_bytes(): swapping bytes for %c (%s)\n",
	      letter(f),fullname(f));
    falcON::swap_bytes(DATA[value(f)], falcON::size(f), NALL);
  }
}
#endif
//
void bodies::block::add_fields(fieldset b) falcON_THROWING
{
  for(fieldbit f; f; ++f)
    if(b.contain(f)) add_field(f);
}
//
void bodies::block::del_fields(fieldset b) falcON_THROWING
{
  for(fieldbit f; f; ++f) 
    if(b.contain(f)) del_field(f);
}
//
void bodies::block::set_fields(fieldset b) falcON_THROWING
{
  for(fieldbit f; f; ++f) 
    if(b.contain(f)) add_field(f);
    else             del_field(f);
}
//
bodies::block::~block() falcON_THROWING
{
  for(fieldbit f; f; ++f)
    del_field(f);
}
//
bodies::block::block(unsigned no,                  // I: our No
		     unsigned na,                  // I: data to allocate
		     unsigned nb,                  // I: # bodies <= na_b
		     unsigned fst,                 // I: first body index
		     bodytype typ,                 // I: hold sph bodies?
		     fieldset bits,                // I: data to allocate
		     bodies  *bods)                // I: pointer to my bodies
  falcON_THROWING
: TYPE       ( typ ),
  NALL       ( na ), 
  NBOD       ( nb ), 
  NO         ( no ),
  FIRST      ( fst ),
  LOCALFIRST ( fst ),
  NEXT       ( 0 ),
  BODS       ( bods )
{
  if(na<nb)
    falcON_THROW("in bodies::block::block(): N_alloc < N_bodies");
  DebugInfo(6,"bodies::block: na=%d, bits=%s, type=%s allowed bits=%s\n",
	    na, word(bits), TYPE.name(), word(bits&TYPE.allows()));
  bits &= TYPE.allows();
  for(fieldbit f; f; ++f)
    set_data_void(f,0);
  try {
    add_fields(bits);
  } catch(exception& E) {
    del_fields(fieldset::all);
    for(fieldbit f; f; ++f)
      del_field(f);
    falcON_RETHROW(E);
  }
}
//
fieldset bodies::block::copy_body(unsigned fr, unsigned to, fieldset _copy)
{
  if(fr>=NALL)
    falcON_THROW("in bodies::block::copy_body(): "
		 "from=%d > NALL=%d\n", fr,NALL);
  if(to>=NALL)
    falcON_THROW("in bodies::block::copy_body(): "
		 "to=%d > NALL=%d\n", to,NALL);
  fieldset copied(fieldset::empty);
  if(fr!=to) {
    for(fieldbit f; f; ++f)
      if(_copy.contain(f) && data_void(f)) {
	memcpy(static_cast<      char*>(data_void(f))+to*falcON::size(f),
	       static_cast<const char*>(data_void(f))+fr*falcON::size(f),
	       falcON::size(f));
	copied |= fieldset(f);
      }
    DebugInfo(8,"bodies::block::copy_body(): copied %s from %d to %d\n",
	      word(copied),fr,to);
  }
  return copied;
}
//
fieldset bodies::block::copy_bodies(const block*that,
				    unsigned    fr,
				    unsigned    to,
				    unsigned    n,
				    fieldset   _copy) falcON_THROWING
{
  if(this == that)
    falcON_THROW("in bodies::block::copy_bodies() from same block");
  if(to+n > NALL)
    falcON_THROW("in bodies::block::copy_bodies(): "
		 "to+n=%d > NALL=%d\n", to+n,NALL);
  if(fr+n > that->NALL)
    falcON_THROW("in bodies::block::copy_bodies(): "
		 "from+n=%d > that->NALL=%d\n", fr+n,that->NALL);
  fieldset copied(fieldset::empty);
  for(fieldbit f; f; ++f)
    if(_copy.contain(f) && data_void(f) && that->data_void(f)) {
      memcpy(static_cast<      char*>(this->data_void(f))+to*falcON::size(f),
	     static_cast<const char*>(that->data_void(f))+fr*falcON::size(f),
	     n*falcON::size(f));
      copied |= fieldset(f);
    }
  return copied;
}
//
inline void bodies::block::skip(unsigned&from,
				flags    copyflag) const falcON_THROWING
{
  if(copyflag)
    for(; from<NBOD && !(flag(from).are_set(copyflag)); ++from ) {}
}
//
// copy up to NALL bodies
// - we copy only bodies of the same type as hold here
// - we only copy bodies whose flag matches last argument
// - if the block copied is finished, we take its NEXT, starting at i=0
// - on return the block pointer is either NULL or together with i they
//   give the first body which was not copied because NALL was exceeded
// - NBOD is set to the bodies copied
fieldset bodies::block::copy(const block*&From,
			     unsigned    &from,
			     fieldset     copydata,
			     flags        copyflag) falcON_THROWING
{
  if( From == this )
    falcON_THROW("in bodies::block::copy(): cannot copy from self");
  NBOD = 0u;
  if(From == 0) return fieldset::empty;
  unsigned _copy;
  unsigned free = NALL;
  fieldset copied;
  if(copyflag && ! has_field(fieldbit::f) )
    falcON_THROW("in bodies::block::copy(): "
		 "copyflag!=0 but flags not supported");
  // skip bodies not to be copied                                               
  From->skip(from,copyflag);
  while(free &&                              // WHILE  we have still space
	From &&                              //   AND  the copied block is valid
	From->TYPE == TYPE &&                //   AND  its type is ours
	from < From->NBOD ) {                //   AND  the index is valid too
    // determine number of bodies to be copied to position from
    if(copyflag) {
      _copy = 0u;
      for(unsigned to=from;
	  to < From->NBOD && From->flag(to).are_set(copyflag) && _copy<free;
	  ++_copy, ++to) {}
    } else
      _copy = min(free, From->NBOD - from);
    // if any body to be copied, copy data, adjust free, NBOD, from, copied
    if(_copy) {
      fieldset c = copy_bodies(From, from, NBOD, _copy, copydata);
      free -= _copy;
      NBOD += _copy;
      from += _copy;
      copied = copied? c : copied & c;
    }
    // skip bodies not to be copied
    From->skip(from,copyflag);
    // end of input block? then take next block
    if(from == From->NBOD) {
      From = From->NEXT;
      if(From == this) 
	falcON_THROW("in bodies::block::copy(): cannot copy from self");
      from = 0;
      if(From) From->skip(from,copyflag);
    }
  }
  return copied;
}
//
void bodies::block::remove(unsigned &removed) falcON_THROWING
{
  if(NBOD == 0) return;
  if(0 == DATA[fieldbit::f] )
    falcON_THROW("in bodies::remove(): flags needed but not supported");
  unsigned lo=0u, hi=NBOD-1;
  while(lo < hi) {
    while(! to_remove(const_datum<fieldbit::f>(lo)) && lo < hi) ++lo;
    while(  to_remove(const_datum<fieldbit::f>(hi)) && lo < hi) --hi;
    if(lo < hi) copy_body(hi--,lo++);
  }
  if(lo == hi && ! to_remove(const_datum<fieldbit::f>(lo))) ++lo;
  removed += NBOD - lo;
  NBOD     = lo;
  DebugInfo(6,"bodies::block::remove(): removed %d: NBOD=%d\n",removed,NBOD);
}
//
// class falcON::bodies
//

// link a new block in
void bodies::add_block(block*B)
{
  // link to last block of same or earlier type, if any, otherwise make it first
  block**P=&FIRST;
  while(*P && (*P)->TYPE <= B->TYPE) P = &((*P)->NEXT);
  B->link(*P);
  *P = B;
  // update TYPES[]
  if(0==TYPES[B->type()])
    TYPES[B->type()]=B;
  // update BLOCK[] and block::NO
  for(int I=0; I!=index::max_blocks; ++I)
    if(BLOCK[I] == 0) {
      BLOCK[I] = B;
      B->NO    = I;
      break;
    }
  // update NBLK and block::BODS
  B->BODS  = this;
  NBLK ++;
  // update block::FIRST and NALL[], NBOD[], NTOT
  set_firsts();
}
// erase a block from our linkage
void bodies::erase_block(block*B)
{
  if(B==0) return;
  // remove from FIRST
  if(FIRST == B)
    FIRST = B->next();
  // remove from TYPES[]
  if(TYPES[B->type()] == B)
    TYPES[B->type()] = B->next_of_same_type();
  // remove from block::NEXT
  for(int i=0; i!=index::max_blocks; ++i)
    if(BLOCK[i] && BLOCK[i]->next() == B) {
      BLOCK[i]->link(B->next());
      break;
    }
  // remove from BLOCK[]
  bool found = false;
  for(int i=0; i!=index::max_blocks; ++i)
    if(BLOCK[i] == B) {
      BLOCK[i] = 0;
      found = true;
      break;
    }
  // reset # bodies info and block::FIRST
  if(found) {
    --NBLK;
    B->BODS = 0;
    set_firsts();
  } else
    falcON_Warning("bodies::erase_block(): block not found in table\n");
}
// remove empty blocks
void bodies::remove_empty_blocks(bool all) falcON_THROWING
{
  for(;;) {
    block*B=0;
    // search for empty block
    for(int i=0; i!=index::max_blocks; ++i)
      if(BLOCK[i] && (all? BLOCK[i]->N_alloc():BLOCK[i]->N_bodies())==0) {
	B=BLOCK[i];
	break;
      }
    // found one: remove it
    if(B) {
      erase_block(B);
      falcON_DEL_O(B);
    } else
      break;
  }
}
// create a new block and link it in
bodies::block* bodies::new_block(bodytype t, unsigned Na, unsigned Nb,
				 fieldset f) falcON_THROWING
{
  if(Nb > Na)
    falcON_THROW("bodies::new_block(): Nb=%u > Na=%u\n",Nb,Na);
  if(Na > index::max_bodies) 
    falcON_THROW("bodies::new_block(): asked for %u > %u bodies\n",
		 Na,  index::max_bodies);
  if(NBLK >= index::max_blocks)
    falcON_THROW("bodies::new_block(): number of blocks exceeded\n");
  block*B=new block(0,Na,Nb,0,t,f,this);
  NNEW[t] += Nb;
  add_block(B);
  DebugInfo(2,"bodies::new_block(): created block for up to %u bodies "
	    "(%u active) of type %s\n", Na,Nb,t.name());
  return B;
}
// reset blocks' FIRST entries (used in parallel code)
void bodies::reset_firsts(unsigned fst[bodytype::NUM])
{
  for(bodytype t; t; ++t) {
    unsigned L=0;
    for(block*B=TYPES[t]; B; B=B->next_of_same_type()) {
      B->set_first(L+fst[t], L);
      L+=B->N_bodies();
    }
  }
}
// # bodies not flagged to be ignored
unsigned bodies::N_subset() const
{
  if(!have(fieldbit::f)) return N_bodies();
  unsigned n = 0;
  LoopAllBodies(this,b) if(in_subset(b)) ++n;
  return n;
}
// delete all blocks and reset related data
void bodies::del_data() falcON_THROWING
{
  for(unsigned i=0; i!=index::max_blocks; ++i) {
    if(BLOCK[i]) falcON_DEL_O(BLOCK[i]);
    BLOCK[i] = 0;
  }
  NBLK = 0u;
  for(bodytype t; t; ++t) {
    NALL [t] = 0u;
    NBOD [t] = 0u;
    TYPES[t] = 0;
  }
  NTOT  = 0;
  FIRST = 0;
}
// destruction: delete all data
bodies::~bodies() falcON_THROWING
{
  DebugInfo(6,"bodies::~bodies(): destructing bodies");
  BITS = fieldset::empty;
  if(C_FORTRAN)
    for(fieldbit f; f; ++f)
      const_cast<block*>(FIRST)->set_data_void(f,0);
  del_data();
}
// set block::FIRST, LOCALFIRST, NALL, NBOD & NTOT
void bodies::set_firsts()
{
  for(bodytype t; t; ++t) {
    NALL[t] = 0u;
    NBOD[t] = 0u;
  }
  NTOT = 0u;
  for(block*P=FIRST; P; P=P->next()) {
    P->set_first(NTOT);
    NALL[P->type()] += P->N_alloc ();
    NBOD[P->type()] += P->N_bodies();
    NTOT            += P->N_bodies();
  }
}
// set up blocks to hold N[t] bodies of type t
void bodies::set_data(const unsigned N[bodytype::NUM]) falcON_THROWING
{
  DebugInfo(5,"bodies::set_data(): N=[%d,%d,%d], BITS=%s\n",
	    N[0],N[1],N[2],word(BITS));
  del_data();
  try {
    block   *last = 0;
    unsigned i    = 0;
    for(bodytype t; t; ++t) {
      NBOD[t]  = NALL[t] = N[t];
      NTOT    += NBOD[t];
      NDEL[t]  = 0u;
      NNEW[t]  = 0u;
      TYPES[t] = 0;
      for(unsigned a,n=0u; n < NALL[t]; n+=a) {
	if(NBLK == index::max_blocks)
	  falcON_THROW("bodies: # blocks exceeds limit");
	a = min(NALL[t]-n, unsigned(index::max_bodies));
	block *b = new block(NBLK,a,a,i,t,BITS,this);
	DebugInfo(10,"allocated %s @ %p\n",
		  nameof(block),static_cast<void*>(b));
	i += a;
	if(last) last->link(b);
	last = b;
	if(n==0u) TYPES[t] = b;
	BLOCK[NBLK++] = b;
      }
    }
  } catch(falcON::exception& E) {
    del_data();
    falcON_RETHROW(E);
  }
  FIRST = BLOCK[0];
  DebugInfo(6,"bodies::set_data(): done\n");
}
// construction 0: construction with N=0, but data fields
bodies::bodies(fieldset bits) falcON_THROWING : 
  BITS      ( bits ),
  C_FORTRAN ( 0 ),
  FORCES    ( 0 )
{
  unsigned n[bodytype::NUM]={0u};
  DebugInfo(2,"bodies::bodies(): constructing bodies @%p: n=%u,%u,%u, bits=%s",
	    this,n[0],n[1],n[2],word(BITS));
  for(unsigned i=0; i!=index::max_blocks; ++i) BLOCK[i] = 0;
  set_data(n);
  set_firsts();
  DebugInfo(2,"bodies::bodies(): constructed\n");
}
// construction 1, new version
bodies::bodies(const unsigned n[bodytype::NUM],
	       fieldset       bits) falcON_THROWING : 
  BITS      ( bits ),
  C_FORTRAN ( 0 ),
  FORCES    ( 0 )
{
  DebugInfo(2,"bodies::bodies(): constructing bodies @%p: n=%u,%u,%u, bits=%s",
	    this,n[0],n[1],n[2],word(bits));
  for(unsigned i=0; i!=index::max_blocks; ++i) BLOCK[i] = 0;
  set_data(n);
  set_firsts();
}
// resets N, data; same as destruction followed by constructor 1            
void bodies::reset(const unsigned n[bodytype::NUM],
		   fieldset       bits) falcON_THROWING
{
  bool keepN = true;
  for(bodytype t; t; ++t) keepN = keepN && NALL[t] == n[t];
  if(keepN) {
    NTOT = 0u;
    for(bodytype t; t; ++t) {
      NBOD[t] = NALL[t];
      NDEL[t] = 0u;
      NNEW[t] = 0u;
      NTOT   += NBOD[t];
    }
    for(unsigned i=0; i!=index::max_blocks; ++i) 
      if(BLOCK[i]) BLOCK[i]->NBOD = BLOCK[i]->NALL;
    del_fields(BITS - bits);
    add_fields(bits - BITS);
  } else {
    del_data();
    BITS = bits;
    set_data(n);
  }
  set_firsts();
}
// construction 2:
// just make a copy of existing bodies:
// - only copy data specified by 2nd argument
// - only copy bodies whose flags matches 3rd argument
// - only copy bodies whose type is contained in 4th argument
bodies::bodies(bodies const&Other,
	       fieldset     copydata,
	       flags        copyflag,
	       bodytypes    copytypes) falcON_THROWING :
  BITS      ( copydata & Other.BITS ),
  C_FORTRAN ( 0 ),
  FORCES    ( 0 )
{
  if(copyflag && !Other.have_flag() ) 
    falcON_THROW("in bodies::bodies(): "
		 "copyflag !=0, but other bodies not supporting flag");
  unsigned n[bodytype::NUM]={0};
  for(bodytype t; t; ++t) if(copytypes.contain(t)) {
    if(copyflag) {
      LoopTypedBodies(&Other,i,t)
	  if( falcON::flag(i).are_set(copyflag) ) ++(n[t]);
    } else 
      n[t] = Other.NBOD[t];
  }
  for(unsigned i=0; i!=index::max_blocks; ++i) BLOCK[i] = 0;
  set_data(n);
  for(bodytype t; t; ++t) if(TYPES[t]) {
    block      *p =TYPES[t];
    const block*op=Other.TYPES[t];
    unsigned    oi=0;
    while(p && op && oi < op->N_bodies()) {
      p->copy(op,oi,copydata,copyflag);
      p = p->next();
    }
  }
  set_firsts();
}
// construction for C & FORTRAN support                                     
bodies::bodies(char, const unsigned n[bodytype::NUM]) falcON_THROWING
: BITS      ( fieldset::empty ),
  C_FORTRAN ( 1 ),
  FORCES    ( 0 )
{
  DebugInfo(3,"bodies::bodies(): constructing bodies for C & FORTRAN: n=%u,%u",
	    n[0],n[1]);
  for(bodytype t; t; ++t)
    if(n[t] > index::max_bodies)
      falcON_THROW("too many bodies\n");
  for(unsigned i=0; i!=index::max_blocks; ++i) BLOCK[i] = 0;
  set_data(n);
  set_firsts();
}
//
void bodies::reset(char, fieldbit f, void*D) falcON_THROWING
{
  if(!C_FORTRAN || !FIRST || BLOCK[0] != FIRST)
    falcON_THROW("bodies::reset() called from wrongly initialized bodies");
  if(D) {
    char* DATA = static_cast<char*>(D);
    BITS |= fieldset(f);
    for(bodytype t; t; ++t)
      if( TYPES[t] && t.allows(f) ) {
	TYPES[t]->set_data_void(f,DATA);
	DATA += falcON::size(f) * TYPES[t]->NALL;
      }
  }
}
//
#if 0
void bodies::swap_bytes(fieldbit f) falcON_THROWING
{
  if(!BITS.contain(f))
    for(const block*p=FIRST; p; p=p->next())
      const_cast<block*>(p)->swap_bytes(f);
}
#endif
//
void bodies::add_field(fieldbit f) falcON_THROWING
{
  if(!BITS.contain(f)) {
    for(const block*p=FIRST; p; p=p->next())
      const_cast<block*>(p)->add_field(f);
    BITS |= fieldset(f);
    if(f == fieldbit::k) reset_keys();
  }
}
//
void bodies::add_fields(fieldset b) falcON_THROWING
{
  if(!BITS.contain(b)) {
    for(const block *p=FIRST; p; p=p->next())
      const_cast<block*>(p)->add_fields(b);
    if(!BITS.contain(fieldbit::k) && b.contain(fieldbit::k)) reset_keys();
    BITS |= b;
  }
}
//
void bodies::del_field(fieldbit f) falcON_THROWING
{
  for(const block *p=FIRST; p; p=p->next())
    const_cast<block*>(p)->del_field(f);
  BITS &= ~(fieldset(f));
}
//
void bodies::del_fields(fieldset b) falcON_THROWING
{
  for(const block *p=FIRST; p; p=p->next())
    const_cast<block*>(p)->del_fields(b);
  BITS &= ~b;
}
//
void bodies::remove(bodytype t)
{
  for(block*P=TYPES[t]; P && P->TYPE == t; P=P->NEXT)
    P->remove(NDEL[t]);
  set_firsts();
  DebugInfo(5,"bodies::remove(%s): removed %d bodies\n", t.name(), NDEL[t]);
}
//
void bodies::remove() {
  for(block*P=FIRST; P; P=P->NEXT)
    P->remove(NDEL[P->TYPE]);
  set_firsts();
  DebugInfo(5,"bodies::remove(): removed %d,%d,%d bodies\n",
	    NDEL[0],NDEL[1],NDEL[2]);
}
//
//
void bodies::merge(bodies&Other) falcON_THROWING
{
  if(NBLK + Other.NBLK > index::max_blocks)
    falcON_THROW("bodies::merge(): too many blocks\n");
  // add blocks
  for(block*P=Other.FIRST; P; P=P->next())
    add_block(P);
  // reset Other (so that its destructor will not harm us)
  Other.FIRST = 0;
  for(bodytype t; t; ++t) {
    Other.TYPES[t] = 0;
    Other.NALL [t] = 0;
    Other.NBOD [t] = 0;
    Other.NNEW [t] = 0;
    Other.NDEL [t] = 0;
  }
  Other.BITS = fieldset::empty;
  Other.NTOT = 0;
  Other.NBLK = 0;
  for(int i=0; i!=index::max_blocks; ++i)
    Other.BLOCK[i] = 0;
}
//
namespace {
  bodies *CopyFrom, *CopyTo;
  Array<bodies::index> IndexTable;
  template<int BIT> struct CopyInOrder {
    static void act(bodytype t)
    {
      unsigned i=0;
      LoopTypedBodies(CopyTo,b,t)
	b.datum<BIT>() = CopyFrom->const_datum<BIT>(IndexTable[i++]);
    }
  };
}
//
//
bodies::block* bodies::ensure_contiguous(unsigned N, bodytype t, unsigned Na)
{
  // do we have a contiguous set of at least N free bodies?
  block*P=TYPES[t];
  while(P && P->N_free()==0) P=P->next_of_same_type();
  unsigned Nf=0;
  for(block*B=P; B; B=B->next_of_same_type()) {
    if     (P==B)       Nf = B->N_free ();       // first part
    else if(B->NBOD==0) Nf+= B->N_alloc();       // further part
    else {                                       // oops: non-contiguous
      while(B && B->N_free()==0) B=B->next_of_same_type();
      Nf= B? B->N_free() : 0;
      P = B;
    }
    if(Nf >= N) break;
  }
  // YES: return block with first part
  if(Nf >= N) {
    DebugInfo(5,"bodies::ensure_contiguous(): found contiguous chunk\n");
    return P;
  }
  // NO:  have to make it
  DebugInfo(5,"bodies::ensure_contiguous(): making new block ...\n");
  return new_block(t,max(Na,N),0,BITS);
}
//
bodies::iterator bodies::new_bodies(unsigned N, bodytype t, unsigned Na)
  falcON_THROWING
{
  // ensure we have enough bodies to activate
  block*P=ensure_contiguous(N,t,Na);
  if(0==P || 0==P->N_free())
    falcON_THROW("bodies::new_bodies(): error in ensure_contiguous()\n");
  unsigned n=N;
  iterator I0(P,P->NBOD);
  // activate bodies
  for(block*B=P; n && B; B=B->next_of_same_type()) {
    unsigned s = min(B->N_free(),n);
    B->NBOD += s;
    n       -= s;
  }
  if(n) falcON_THROW("bodies::new_bodies(): cannot find enough free bodies\n");
  set_firsts();
  // flag as new
  if(have(fieldbit::f)) {
    iterator IN(I0,N);
    for(iterator I(I0); I!=IN; ++I)
      I.flag().add(flags::newbody);
  }
  return I0;
}
//
bodies::iterator bodies::new_body(bodytype t, unsigned Na) falcON_THROWING
{
  // ensure we have a body to activate
  block*P=ensure_contiguous(1,t,Na);
  if(0==P || 0==P->N_free())
    falcON_THROW("bodies::new_body(): error in ensure_contiguous()\n");
  // activate body
  iterator I0(P,P->NBOD++);
  set_firsts();
  // flag as new
  if(have(fieldbit::f)) I0.flag().add(flags::newbody);
  return I0;
}
//
void bodies::joinup(bodytype t) falcON_THROWING
{
  block*To=TYPES[t];
  bool copy=false;
  for(;;) {
    // find block of type t with free space
    while(To && To->N_free()==0) To=To->next_of_same_type();
    if(0 == To) break;
    // find later block of type t with data
    block*Fr=To->next_of_same_type();
    while(Fr && Fr->NBOD==0) Fr=Fr->next_of_same_type();
    if(0 == Fr) break;
    // copy data from Fr to To
    int nc=min(To->N_free(), Fr->NBOD);
    To->copy_bodies(Fr,Fr->NBOD-nc,To->NBOD,nc);
    To->NBOD += nc;
    Fr->NBOD -= nc;
    copy = true;
  }
  if(copy) set_firsts();
}
//
falcON::real bodies::TotalMass(bodytype t) const
{
  if(!t || TYPES[t]==0 || !(TYPES[t]->has_field(fieldbit::m)) )
    return zero;
  real M(zero);
  for(const block* b=TYPES[t]; b; b=b->next_of_same_type())
    for(unsigned i=0; i!=b->N_bodies(); ++i)
      M += b->const_datum<fieldbit::m>(i);
  return M;
}
// sorted index table
void bodies::sorted(Array<index>&table, 
		    real       (*func)(iterator const&)) const falcON_THROWING
{
  const int n = N_subset();
  real *Q = falcON_NEW(real, n);
  index*I = falcON_NEW(index,n);
  if(have(fieldbit::f)) {
    int i = 0;
    LoopSubsetBodies(this,b) {
      I[i] = index(b);
      Q[i] = func(b);
      ++i;
    }
  } else {
    int i = 0;
    LoopAllBodies(this,b) {
      I[i] = index(b);
      Q[i] = func(b);
      ++i;
    }
  }
  int*R = falcON_NEW(int,n);
  HeapIndex(Q,n,R);
  table.reset(n);
  for(int i=0; i!=n; ++i)
    table[i] = I[R[i]];
  falcON_DEL_A(Q);
  falcON_DEL_A(I);
  falcON_DEL_A(R);
}
//
void bodies::sorted(Array<index>&table, 
		    Array<real> &quant, 
		    real       (*func)(iterator const&)) const falcON_THROWING
{
  const int n = N_subset();
  real *Q = falcON_NEW(real, n);
  index*I = falcON_NEW(index,n);
  if(have(fieldbit::f)) {
    int i = 0;
    LoopSubsetBodies(this,b) {
      I[i] = index(b);
      Q[i] = func(b);
      ++i;
    }
  } else {
    int i = 0;
    LoopAllBodies(this,b) {
      I[i] = index(b);
      Q[i] = func(b);
      ++i;
    }
  }
  int*R = falcON_NEW(int,n);
  HeapIndex(Q,n,R);
  table.reset(n);
  quant.reset(n);
  for(int i=0; i!=n; ++i) {
    table[i] = I[R[i]];
    quant[i] = Q[R[i]];
  }
  falcON_DEL_A(Q);
  falcON_DEL_A(I);
  falcON_DEL_A(R);
}
//
namespace {
  struct Nbour { real Q; bodies::index I; };
  inline bool operator<(Nbour const&a, Nbour const&b) {return a.Q<b.Q;}
//   inline bool operator>(Nbour const&a, Nbour const&b) {return a.Q>b.Q;}
//   inline bool operator<(real q, Nbour const&b) {return q<b.Q;}
//   inline bool operator>(real q, Nbour const&b) {return q>b.Q;}
//   inline bool operator<(Nbour const&a, real q) {return a.Q<q;}
//   inline bool operator>(Nbour const&a, real q) {return a.Q>q;}
  real Huge = 1.e30;
}
falcON_TRAITS(Nbour,"<anonymous>::Nbour");
//
unsigned bodies::findNeighbours(const body&B, unsigned K, Array<index>&I) const
  falcON_THROWING
{ 
  if(!have_pos())
    falcON_THROW("bodies::findNeighbours(): have no positions\n");
  Nbour*List = falcON_NEW(Nbour,K);
  for(unsigned k=0; k!=K; ++k)
    List[k].Q = Huge;
  unsigned Niac=0;
  LoopSubsetBodies(this,b) {
    real q = dist_sq(falcON::pos(B),falcON::pos(b));
    if(List->Q > q) {
      List->Q = q;
      List->I = bodies::index(b);
      MaxHeap::after_top_replace(List,K);
      ++Niac;
    }
  }
  MaxHeap::sort(List,K);
  I.reset(K);
  if(Niac < K) K = Niac;
  for(unsigned k=0; k!=K; ++k)
    I[k] = List[k].I;
  falcON_DEL_A(List);
  return K;
}
