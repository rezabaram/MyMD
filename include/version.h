// This file is a part of Molecular Dynamics code for
// simulating ellipsoidal packing. The author cannot
// guarantee the correctness nor the intended functionality.
//
// March 2012, Reza Baram


#ifndef MYMD_VERSION_H
#define MYMD_VERSION_H

// CMake defines this from project(MyMD VERSION ...), so the two cannot drift.
// The fallback is for the Makefile build, which has no such notion.
#ifndef MYMD_VERSION
#define MYMD_VERSION "0.1.0"
#endif

#endif /* MYMD_VERSION_H */
