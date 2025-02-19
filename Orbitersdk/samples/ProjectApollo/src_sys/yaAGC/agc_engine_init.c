/*
  Copyright 2003-2005 Ronald S. Burkey <info@sandroid.org>

  This file is part of yaAGC.

  yaAGC is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 2 of the License, or
  (at your option) any later version.

  yaAGC is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with yaAGC; if not, write to the Free Software
  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

  In addition, as a special exception, Ronald S. Burkey gives permission to
  link the code of this program with the Orbiter SDK library (or with 
  modified versions of the Orbiter SDK library that use the same license as 
  the Orbiter SDK library), and distribute linked combinations including 
  the two. You must obey the GNU General Public License in all respects for 
  all of the code used other than the Orbiter SDK library. If you modify 
  this file, you may extend this exception to your version of the file, 
  but you are not obligated to do so. If you do not wish to do so, delete 
  this exception statement from your version. 
 
  Filename:	agc_engine_init.c
  Purpose:	This is the function which initializes the AGC simulation,
  		from a file representing the binary image of core memory.
  Compiler:	GNU gcc.
  Contact:	Ron Burkey <info@sandroid.org>
  Reference:	http://www.ibiblio.org/apollo/index.html
  Mods:		04/05/03 RSB.	Began.
  		09/07/03 RSB.	Fixed data ordering in the core-rope image
				file (to both endian CPU types to work).
		11/26/03 RSB.	Up to now, a pseudo-linear space was used to
				model internal AGC memory.  This was simply too
				tricky to work with, because it was too hard to
				understand the address conversions that were
				taking place.  I now use a banked model much
				closer to the true AGC memory map.
		11/29/03 RSB.	Added the core-dump save/load.
		05/06/04 RSB	Now use rfopen in looking for the binary.
		07/12/04 RSB	Q is now 16 bits.
		07/15/04 RSB	AGC data now aligned at bit 0 rathern then 1.
		07/17/04 RSB	I/O channels 030-033 now default to 077777
				instead of 00000, since the signals are 
				supposed to be inverted.
		02/27/05 RSB	Added the license exception, as required by
				the GPL, for linking to Orbiter SDK libraries.
		05/14/05 RSB	Corrected website references.
	 	07/05/05 RSB	Added AllOrErasable.
		07/07/05 RSB	On a resume, now restores 010 on up (rather
				than 020 on up), on Hugh's advice.
		09/30/16 MAS	Added initialization of NightWatchman.
		01/04/17 MAS    Added initialization of ParityFail.
		01/30/17 MAS    Added support for heuristic loading of ROM files
						produced with --hardware, by looking for any set
		                parity bits. If such a file is detected, parity
		                bit checking is enabled.
		03/09/17 MAS    Added initialization of SbyStillPressed.
		03/26/17 MAS    Added initialization of previously-static things
		                 from agc_engine.c that are now in agc_t.
		03/27/17 MAS    Fixed a parity-related program loading bug and
                         added initialization of a new night watchman bit.
		04/02/17 MAS	Added initialization of a couple of flags used
						for simulation of the TC Trap hardware bug.
		04/16/17 MAS    Added initialization of warning filter variables.
		05/16/17 MAS    Enabled interrupts at startup.
		07/13/17 MAS	Added initialization of the three HANDRUPT traps.
		01/29/24 MAS	Added initialization of RadarGateCounter.
*/

#include <stdio.h>
#include <string.h>
#include "agc_engine.h"
FILE *rfopen (const char *Filename, const char *mode);

//---------------------------------------------------------------------------
// Returns:
//      0 -- success.
//      1 -- ROM image file not found.
//      2 -- ROM image file larger than core memory.
//      3 -- ROM image file size is odd.
//      4 -- agc_t structure not allocated.
//      5 -- File-read error.
//      6 -- Core-dump file not found.
// Normally, on input the CoreDump filename is NULL, in which case all of the 
// i/o channels, erasable memory, etc., are cleared to their reset values.
// When the CoreDump is loaded instead, it allows execution to continue precisely
// from the point at which the CoreDump was created, if AllOrErasable != 0.
// If AllOrErasable == 0, then only the erasable memory is initialized from the
// core-dump file.

int
agc_load_binfile(agc_t *State, const char *RomImage)

{
  FILE *fp = NULL;
  int Bank;
  int m, n, i, j;

  // The following sequence of steps loads the ROM image into the simulated
  // core memory, in what I think is a pretty obvious way.

  int RetVal = 1;
  fp = rfopen (RomImage, "rb");
  if (fp == NULL)
    goto Done;

  RetVal = 3;
  fseek (fp, 0, SEEK_END);
  n = ftell (fp);
  if (0 != (n & 1))		// Must be an integral number of words.
    goto Done;

  RetVal = 2;
  n /= 2;			// Convert byte-count to word-count.
  if (n > 36 * 02000)
    goto Done;

  RetVal = 4;
  fseek (fp, 0, SEEK_SET);
  if (State == NULL)
    goto Done;

  State->CheckParity = 0;
  memset(&State->Parities, 0, sizeof(State->Parities));

  RetVal = 5;
  Bank = 2;
  for (Bank = 2, j = 0, i = 0; i < n; i++)
    {
      unsigned char In[2];
	  uint8_t Parity;
	  uint16_t RawValue;
      m = fread (In, 1, 2, fp);
      if (m != 2)
	goto Done;
      // Within the input file, the fixed-memory banks are arranged in the order
      // 2, 3, 0, 1, 4, 5, 6, 7, ..., 35.  Therefore, we have to take a little care
      // reordering the banks.
      if (Bank > 35)
	{
	  RetVal = 2;
	  goto Done;
	}
	  RawValue = (In[0] * 256 + In[1]);
	  Parity = RawValue & 1;

	  State->Fixed[Bank][j] = RawValue >> 1;
	  State->Parities[(Bank * 02000 + j) / 32] |= Parity << (j % 32);
	  j++;

	  // If any of the parity bits are actually set, this must be a ROM built with
	  // --hardware. Enable parity checking.
	  if (Parity)
	    State->CheckParity = 1;

      if (j == 02000)
	{
	  j = 0;
	  // Bank filled.  Advance to next fixed-memory bank.
	  if (Bank == 2)
	    Bank = 3;
	  else if (Bank == 3)
	    Bank = 0;
	  else if (Bank == 0)
	    Bank = 1;
	  else if (Bank == 1)
	    Bank = 4;
	  else
	    Bank++;
	}
    }

Done:
  if (fp != NULL)
    fclose (fp);
  return (RetVal);
}

int
agc_engine_init (agc_t * State, const char *RomImage, const char *CoreDump,
		 int AllOrErasable)
{
  int RetVal = 0, i, j, Bank;
  FILE *cd = NULL;

#ifndef WIN32
  // The purpose of this is to make sure that getchar doesn't halt the program
  // when there's no keystroke immediately available.
  UnblockSocket (fileno (stdin));
#endif

  if (RomImage)
	  RetVal = agc_load_binfile(State, RomImage);
 
  // Clear i/o channels.
  for (i = 0; i < NUM_CHANNELS; i++)
    State->InputChannel[i] = 0;
  State->InputChannel[030] = 037777;
  State->InputChannel[031] = 077777;
  State->InputChannel[032] = 077777;
  State->InputChannel[033] = 077777;

  // Clear erasable memory.
  for (Bank = 0; Bank < 8; Bank++)
    for (j = 0; j < 0400; j++)
      State->Erasable[Bank][j] = 0;
  State->Erasable[0][RegZ] = 04000;	// Initial program counter.

  // Set up the CPU state variables that aren't part of normal memory.
  RetVal = 0;
  State->CycleCounter = 0;
  State->ExtraCode = 0;
  State->AllowInterrupt = 1; // The GOJAM sequence enables interrupts
  State->InterruptRequests[8] = 1;	// DOWNRUPT.
  //State->RegA16 = 0;
  State->PendFlag = 0;
  State->PendDelay = 0;
  State->ExtraDelay = 0;
  //State->RegQ16 = 0;

  State->NightWatchman = 0;
  State->NightWatchmanTripped = 0;
  State->RuptLock = 0;
  State->NoRupt = 0;
  State->TCTrap = 0;
  State->NoTC = 0;
  State->ParityFail = 0;

  State->WarningFilter = 0;
  State->GeneratedWarning = 0;

  State->RestartLight = 0;
  State->Standby = 0;
  State->SbyPressed = 0;
  State->SbyStillPressed = 0;

  State->NextZ = 0;
  State->ScalerCounter = 0;
  State->ChannelRoutineCount = 0;

  State->DskyTimer = 0;
  State->DskyFlash = 0;
  State->DskyChannel163 = 0;

  State->TookBZF = 0;
  State->TookBZMF = 0;

  State->Trap31A = 0;
  State->Trap31B = 0;
  State->Trap32 = 0;

  State->RadarGateCounter = 0;

  if (CoreDump != NULL)
    {
      cd = fopen (CoreDump, "r");
      if (cd == NULL)
        {
	  if (AllOrErasable)
	    RetVal = 6;
	  else
	    RetVal = 0;
	}
      else
	{
	  RetVal = 5;

	  // Load up the i/o channels.
	  for (i = 0; i < NUM_CHANNELS; i++)
	    {
	      if (1 != fscanf (cd, "%o", &j))
		goto Done;
	      if (AllOrErasable)
	        State->InputChannel[i] = j;
	    }

	  // Load up erasable memory.
	  for (Bank = 0; Bank < 8; Bank++)
	    for (j = 0; j < 0400; j++)
	      {
		if (1 != fscanf (cd, "%o", &i))
		  goto Done;
		if (AllOrErasable || Bank > 0 || j >= 010)
		  State->Erasable[Bank][j] = i;
	      }

          if (AllOrErasable)
	    {
	      // Set up the CPU state variables that aren't part of normal memory.
	      if (1 != fscanf (cd, "%o", &i))
		goto Done;
	      State->CycleCounter = i;
	      if (1 != fscanf (cd, "%o", &i))
		goto Done;
	      State->ExtraCode = i;
	      // I've seen no indication so far of a reset value for interrupt-enable. 
	      if (1 != fscanf (cd, "%o", &i))
		goto Done;
	      State->AllowInterrupt = i;
	      //if (1 != fscanf (cd, "%o", &i))
	      //  goto Done;
	      //State->RegA16 = i;
	      if (1 != fscanf (cd, "%o", &i))
		goto Done;
	      State->PendFlag = i;
	      if (1 != fscanf (cd, "%o", &i))
		goto Done;
	      State->PendDelay = i;
	      if (1 != fscanf (cd, "%o", &i))
		goto Done;
	      State->ExtraDelay = i;
	      //if (1 != fscanf (cd, "%o", &i))
	      //  goto Done;
	      //State->RegQ16 = i;
	    }

	  RetVal = 0;
	}
    }

Done:
  if (cd != NULL)
    fclose (cd);
  return (RetVal);
}

//-------------------------------------------------------------------------------
// A function for creating a core-dump which can be read by agc_engine_init.

void
MakeCoreDump (agc_t * State, const char *CoreDump)
{
#if defined (WIN32) || defined (__APPLE__)  
  uint64_t lli;
#else
  unsigned long long lli;
#endif
  FILE *cd = NULL;
  int i, j, Bank;

  cd = fopen (CoreDump, "w");
  if (cd == NULL)
    {
      printf ("Could not create the core-dump file.\n");
      return;
    }

  // Write out the i/o channels.
  for (i = 0; i < NUM_CHANNELS; i++)
    fprintf (cd, "%06o\n", State->InputChannel[i]);

  // Write out the erasable memory.
  for (Bank = 0; Bank < 8; Bank++)
    for (j = 0; j < 0400; j++)
      fprintf (cd, "%06o\n", State->Erasable[Bank][j]);

  // Write out CPU state variables that aren't part of normal memory.
  fprintf (cd, "%llo\n", State->CycleCounter);
  fprintf (cd, "%o\n", State->ExtraCode);
  fprintf (cd, "%o\n", State->AllowInterrupt);
  //fprintf (cd, "%o\n", State->RegA16);
  fprintf (cd, "%o\n", State->PendFlag);
  fprintf (cd, "%o\n", State->PendDelay);
  fprintf (cd, "%o\n", State->ExtraDelay);
  //fprintf (cd, "%o\n", State->RegQ16);

  // 03/27/09 RSB.  Extra agc_t fields that weren't being saved before
  // now.  I've made no analysis to determine that all of these are 
  // actually needed for anything.
  fprintf  (cd, "%o\n", State->OutputChannel7);
  for (i = 0; i < 16; i++)
    fprintf (cd, "%o\n", State->OutputChannel10[i]);
  fprintf (cd, "%o\n", State->IndexValue);
  for (i = 0; i < 1 + NUM_INTERRUPT_TYPES; i++)
    fprintf (cd, "%o\n", State->InterruptRequests[i]);
  fprintf (cd, "%o\n", State->InIsr);
  fprintf (cd, "%o\n", State->SubstituteInstruction);
  fprintf (cd, "%o\n", State->DownruptTimeValid);
  fprintf (cd, "%llo\n", lli = State->DownruptTime);
  fprintf (cd, "%o\n", 0);

  printf ("Core-dump file \"%s\" created.\n", CoreDump);
  fclose (cd);
  return;

//Done:
//  printf ("Write-error on core-dump file.\n");
//  fclose (cd);
//  return;

}
