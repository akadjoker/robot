/*************************************************************************
 *                                                                       *
 * Open Dynamics Engine, Copyright (C) 2001,2002 Russell L. Smith.       *
 * All rights reserved.  Email: russ@q12.org   Web: www.q12.org          *
 *                                                                       *
 * This library is free software; you can redistribute it and/or         *
 * modify it under the terms of EITHER:                                  *
 *   (1) The GNU Lesser General Public License as published by the Free  *
 *       Software Foundation; either version 2.1 of the License, or (at  *
 *       your option) any later version. The text of the GNU Lesser      *
 *       General Public License is included with this library in the     *
 *       file LICENSE.TXT.                                               *
 *   (2) The BSD-style license that is included with this library in     *
 *       the file LICENSE-BSD.TXT.                                       *
 *                                                                       *
 * This library is distributed in the hope that it will be useful,       *
 * but WITHOUT ANY WARRANTY; without even the implied warranty of        *
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the files    *
 * LICENSE.TXT and LICENSE-BSD.TXT for more details.                     *
 *                                                                       *
 *************************************************************************/

#ifndef _ODE_UTIL_H_
#define _ODE_UTIL_H_

#include "objects.h"
#include "common.h"


/* utility */

void dInternalHandleAutoDisabling (dxWorld *world, dReal stepsize);
void dxStepBody (dxBody *b, dReal h);


struct dxWorldProcessMemoryManager:
    public dBase
{
    typedef void *(*alloc_block_fn_t)(sizeint block_size);
    typedef void *(*shrink_block_fn_t)(void *block_pointer, sizeint block_current_size, sizeint block_smaller_size);
    typedef void (*free_block_fn_t)(void *block_pointer, sizeint block_current_size);

    dxWorldProcessMemoryManager(alloc_block_fn_t fnAlloc, shrink_block_fn_t fnShrink, free_block_fn_t fnFree)
    {
        Assign(fnAlloc, fnShrink, fnFree);
    }

    void Assign(alloc_block_fn_t fnAlloc, shrink_block_fn_t fnShrink, free_block_fn_t fnFree)
    {
        m_fnAlloc = fnAlloc;
        m_fnShrink = fnShrink;
        m_fnFree = fnFree;
    }

    alloc_block_fn_t m_fnAlloc;
    shrink_block_fn_t m_fnShrink;
    free_block_fn_t m_fnFree;
};

extern dxWorldProcessMemoryManager g_WorldProcessMallocMemoryManager;

struct dxWorldProcessMemoryReserveInfo:
    public dBase
{
    dxWorldProcessMemoryReserveInfo(float fReserveFactor, unsigned uiReserveMinimum)
    {
        Assign(fReserveFactor, uiReserveMinimum);
    }

    void Assign(float fReserveFactor, unsigned uiReserveMinimum)
    {
        m_fReserveFactor = fReserveFactor;
        m_uiReserveMinimum = uiReserveMinimum;
    }

    float m_fReserveFactor; // Use float as precision does not matter here
    unsigned m_uiReserveMinimum;
};

extern dxWorldProcessMemoryReserveInfo g_WorldProcessDefaultReserveInfo;


class dxWorldProcessMemArena:
    private dBase // new/delete must not be called for this class
{
public:
#define BUFFER_TO_ARENA_EXTRA (EFFICIENT_ALIGNMENT + dEFFICIENT_SIZE(sizeof(dxWorldProcessMemArena)))
    static bool IsArenaPossible(sizeint nBufferSize)
    {
        return SIZE_MAX - BUFFER_TO_ARENA_EXTRA >= nBufferSize; // This ensures there will be no overflow
    }

    static sizeint MakeBufferSize(sizeint nArenaSize)
    {
        return nArenaSize - BUFFER_TO_ARENA_EXTRA;
    }

    static sizeint MakeArenaSize(sizeint nBufferSize)
    {
        return BUFFER_TO_ARENA_EXTRA + nBufferSize;
    }
#undef BUFFER_TO_ARENA_EXTRA

    bool IsStructureValid() const
    {
        return m_pAllocBegin != NULL && m_pAllocEnd != NULL && m_pAllocBegin <= m_pAllocEnd 
            && (m_pAllocCurrentOrNextArena == NULL || m_pAllocCurrentOrNextArena == m_pAllocBegin) 
            && m_pArenaBegin != NULL && m_pArenaBegin <= m_pAllocBegin; 
    }

    sizeint GetMemorySize() const
    {
        return (sizeint)m_pAllocEnd - (sizeint)m_pAllocBegin;
    }

    void *SaveState() const
    {
        return m_pAllocCurrentOrNextArena;
    }

    void RestoreState(void *state)
    {
        m_pAllocCurrentOrNextArena = state;
    }

    void ResetState()
    {
        m_pAllocCurrentOrNextArena = m_pAllocBegin;
    }

    void *PeekBufferRemainder() const
    {
        return m_pAllocCurrentOrNextArena;
    }

    void *AllocateBlock(sizeint size)
    {
        void *arena = m_pAllocCurrentOrNextArena;
        m_pAllocCurrentOrNextArena = dOFFSET_EFFICIENTLY(arena, size);
        dIASSERT(m_pAllocCurrentOrNextArena <= m_pAllocEnd);
        dIASSERT(dEFFICIENT_PTR(arena) == arena);
        
        return arena;
    }

    void *AllocateOveralignedBlock(sizeint size, unsigned alignment)
    {
        void *arena = m_pAllocCurrentOrNextArena;
        m_pAllocCurrentOrNextArena = dOFFSET_OVERALIGNEDLY(arena, size, alignment);
        dIASSERT(m_pAllocCurrentOrNextArena <= m_pAllocEnd);

        void *block = dOVERALIGNED_PTR(arena, alignment);
        return block;
    }

    template<typename ElementType>
    ElementType *AllocateArray(sizeint count)
    {
        return (ElementType *)AllocateBlock(count * sizeof(ElementType));
    }

    template<typename ElementType>
    ElementType *AllocateOveralignedArray(sizeint count, unsigned alignment)
    {
        return (ElementType *)AllocateOveralignedBlock(count * sizeof(ElementType), alignment);
    }

    template<typename ElementType>
    void ShrinkArray(ElementType *arr, sizeint oldcount, sizeint newcount)
    {
        dIASSERT(newcount <= oldcount);
        dIASSERT(dOFFSET_EFFICIENTLY(arr, oldcount * sizeof(ElementType)) == m_pAllocCurrentOrNextArena);
        m_pAllocCurrentOrNextArena = dOFFSET_EFFICIENTLY(arr, newcount * sizeof(ElementType));
    }

public:
    static dxWorldProcessMemArena *ReallocateMemArena (
        dxWorldProcessMemArena *oldarena, sizeint memreq, 
        const dxWorldProcessMemoryManager *memmgr, float rsrvfactor, unsigned rsrvminimum);
    static void FreeMemArena (dxWorldProcessMemArena *arena);

    dxWorldProcessMemArena *GetNextMemArena() const { return (dxWorldProcessMemArena *)m_pAllocCurrentOrNextArena; }
    void SetNextMemArena(dxWorldProcessMemArena *pArenaInstance) { m_pAllocCurrentOrNextArena = pArenaInstance; }

private:
    static sizeint AdjustArenaSizeForReserveRequirements(sizeint arenareq, float rsrvfactor, unsigned rsrvminimum);

private:
    void *m_pAllocCurrentOrNextArena;
    void *m_pAllocBegin;
    void *m_pAllocEnd;
    void *m_pArenaBegin;

    const dxWorldProcessMemoryManager *m_pArenaMemMgr;
};

class dxWorldProcessContext:
    public dBase
{
public:
    dxWorldProcessContext();
    ~dxWorldProcessContext();

    void CleanupWorldReferences(dxWorld *pswWorldInstance);

public:
    bool EnsureStepperSyncObjectsAreAllocated(dxWorld *pswWorldInstance);
    dCallWaitID GetIslandsSteppingWait() const { return m_pcwIslandsSteppingWait; }

public:
    dxWorldProcessMemArena *ObtainStepperMemArena();
    void ReturnStepperMemArena(dxWorldProcessMemArena *pmaArenaInstance);

    dxWorldProcessMemArena *ReallocateIslandsMemArena(sizeint nMemoryRequirement, 
        const dxWorldProcessMemoryManager *pmmMemortManager, float fReserveFactor, unsigned uiReserveMinimum);
    bool ReallocateStepperMemArenas(dxWorld *world, unsigned nIslandThreadsCount, sizeint nMemoryRequirement, 
        const dxWorldProcessMemoryManager *pmmMemortManager, float fReserveFactor, unsigned uiReserveMinimum);

private:
    static void FreeArenasList(dxWorldProcessMemArena *pmaExistingArenas);

private:
    void SetIslandsMemArena(dxWorldProcessMemArena *pmaInstance) { m_pmaIslandsArena = pmaInstance; }
    dxWorldProcessMemArena *GetIslandsMemArena() const { return m_pmaIslandsArena; }

    void SetStepperArenasList(dxWorldProcessMemArena *pmaInstance) { m_pmaStepperArenas = pmaInstance; }
    dxWorldProcessMemArena *GetStepperArenasList() const { return m_pmaStepperArenas; }

    inline dxWorldProcessMemArena *GetStepperArenasHead() const;
    inline bool TryExtractingStepperArenasHead(dxWorldProcessMemArena *pmaHeadInstance);
    inline bool TryInsertingStepperArenasHead(dxWorldProcessMemArena *pmaArenaInstance, dxWorldProcessMemArena *pmaExistingHead);

public:
    void LockForAddLimotSerialization();
    void UnlockForAddLimotSerialization();
    void LockForStepbodySerialization();
    void UnlockForStepbodySerialization();

private:
    enum dxProcessContextMutex
    {
        dxPCM_STEPPER_ARENA_OBTAIN,
        dxPCM_STEPPER_ADDLIMOT_SERIALIZE,
        dxPCM_STEPPER_STEPBODY_SERIALIZE,

        dxPCM__MAX
    };

    static const char *const m_aszContextMutexNames[dxPCM__MAX];

private:
    dxWorldProcessMemArena  *m_pmaIslandsArena;
    dxWorldProcessMemArena  *volatile m_pmaStepperArenas;
    dxWorld                 *m_pswObjectsAllocWorld;
    dMutexGroupID           m_pmgStepperMutexGroup;
    dCallWaitID             m_pcwIslandsSteppingWait;
};

struct dxWorldProcessIslandsInfo
{
    void AssignInfo(sizeint islandcount, unsigned int const *islandsizes, dxBody *const *bodies, dxJoint *const *joints)
    {
        m_IslandCount = islandcount;
        m_pIslandSizes = islandsizes;
        m_pBodies = bodies;
        m_pJoints = joints;
    }

    sizeint GetIslandsCount() const { return m_IslandCount; }
    unsigned int const *GetIslandSizes() const { return m_pIslandSizes; }
    dxBody *const *GetBodiesArray() const { return m_pBodies; }
    dxJoint *const *GetJointsArray() const { return m_pJoints; }

private:
    sizeint                  m_IslandCount;
    unsigned int const      *m_pIslandSizes;
    dxBody *const           *m_pBodies;
    dxJoint *const          *m_pJoints;
};

struct dxStepperProcessingCallContext
{
    dxStepperProcessingCallContext(dxWorld *world, dReal stepSize, unsigned stepperAllowedThreads, unsigned lcpAllowedThreads,
        dxWorldProcessMemArena *stepperArena, dxBody *const *islandBodiesStart, dxJoint *const *islandJointsStart): 
        m_world(world), m_stepSize(stepSize), m_stepperArena(stepperArena), m_finalReleasee(NULL), 
        m_islandBodiesStart(islandBodiesStart), m_islandJointsStart(islandJointsStart), m_islandBodiesCount(0), m_islandJointsCount(0),
        m_stepperAllowedThreads(stepperAllowedThreads),
        m_lcpAllowedThreads(lcpAllowedThreads)
    {
    }

    void AssignIslandSelection(dxBody *const *islandBodiesStart, dxJoint *const *islandJointsStart, 
        unsigned islandBodiesCount, unsigned islandJointsCount)
    {
        m_islandBodiesStart = islandBodiesStart;
        m_islandJointsStart = islandJointsStart;
        m_islandBodiesCount = islandBodiesCount;
        m_islandJointsCount = islandJointsCount;
    }

    dxBody *const *GetSelectedIslandBodiesEnd() const { return m_islandBodiesStart + m_islandBodiesCount; }
    dxJoint *const *GetSelectedIslandJointsEnd() const { return m_islandJointsStart + m_islandJointsCount; }

    void AssignStepperCallFinalReleasee(dCallReleaseeID finalReleasee)
    {
        m_finalReleasee = finalReleasee;
    }

    dxWorld                 *const m_world;
    dReal                   const m_stepSize;
    dxWorldProcessMemArena  *m_stepperArena;
    dCallReleaseeID         m_finalReleasee;
    dxBody *const           *m_islandBodiesStart;
    dxJoint *const          *m_islandJointsStart;
    unsigned                m_islandBodiesCount;
    unsigned                m_islandJointsCount;
    unsigned                m_stepperAllowedThreads;
    unsigned                m_lcpAllowedThreads;
};

#define BEGIN_STATE_SAVE(memarena, state) void *state = memarena->SaveState();
#define END_STATE_SAVE(memarena, state) memarena->RestoreState(state)

typedef void (*dstepper_fn_t) (const dxStepperProcessingCallContext *callContext);
typedef unsigned (*dmaxcallcountestimate_fn_t) (unsigned activeThreadCount, unsigned steppingAllowedThreadCount, unsigned lcpAllowedThreadCount);

bool dxProcessIslands (dxWorld *world, const dxWorldProcessIslandsInfo &islandsInfo, 
                       dReal stepSize, dstepper_fn_t stepper, dmaxcallcountestimate_fn_t maxCallCountEstimator);


typedef sizeint (*dmemestimate_fn_t) (dxBody * const *body, unsigned int nb, 
                                     dxJoint * const *_joint, unsigned int _nj);

bool dxReallocateWorldProcessContext (dxWorld *world, dxWorldProcessIslandsInfo &islandsinfo, 
                                      dReal stepsize, dmemestimate_fn_t stepperestimate);

dxWorldProcessMemArena *dxAllocateTemporaryWorldProcessMemArena(
    sizeint memreq, const dxWorldProcessMemoryManager *memmgr/*=NULL*/, const dxWorldProcessMemoryReserveInfo *reserveinfo/*=NULL*/);
void dxFreeTemporaryWorldProcessMemArena(dxWorldProcessMemArena *arena);


template<class ClassType>
inline ClassType *AllocateOnDemand(ClassType *&pctStorage)
{
    ClassType *pctCurrentInstance = pctStorage;

    if (!pctCurrentInstance)
    {
        pctCurrentInstance = new ClassType();
        pctStorage = pctCurrentInstance;
    }

    return pctCurrentInstance;
}


// World stepping working memory object
class dxStepWorkingMemory:
    public dBase
{
public:
    dxStepWorkingMemory(): m_uiRefCount(1), m_ppcProcessingContext(NULL), m_priReserveInfo(NULL), m_pmmMemoryManager(NULL) {}

private:
    friend struct dBase; // To avoid GCC warning regarding private destructor
    ~dxStepWorkingMemory() // Use Release() instead
    {
        delete m_ppcProcessingContext;
        delete m_priReserveInfo;
        delete m_pmmMemoryManager;
    }

public:
    void Addref()
    {
        dIASSERT(~m_uiRefCount != 0);
        ++m_uiRefCount;
    }

    void Release()
    {
        dIASSERT(m_uiRefCount != 0);
        if (--m_uiRefCount == 0)
        {
            delete this;
        }
    }

public:
    void CleanupMemory()
    {
        delete m_ppcProcessingContext;
        m_ppcProcessingContext = NULL;
    }

    void CleanupWorldReferences(dxWorld *world)
    {
        if (m_ppcProcessingContext != NULL)
        {
            m_ppcProcessingContext->CleanupWorldReferences(world);
        }
    }

public: 
    dxWorldProcessContext *SureGetWorldProcessingContext() { return AllocateOnDemand(m_ppcProcessingContext); }
    dxWorldProcessContext *GetWorldProcessingContext() const { return m_ppcProcessingContext; }

    const dxWorldProcessMemoryReserveInfo *GetMemoryReserveInfo() const { return m_priReserveInfo; }
    const dxWorldProcessMemoryReserveInfo *SureGetMemoryReserveInfo() const { return m_priReserveInfo ? m_priReserveInfo : &g_WorldProcessDefaultReserveInfo; }
    void SetMemoryReserveInfo(float fReserveFactor, unsigned uiReserveMinimum)
    {
        if (m_priReserveInfo) { m_priReserveInfo->Assign(fReserveFactor, uiReserveMinimum); }
        else { m_priReserveInfo = new dxWorldProcessMemoryReserveInfo(fReserveFactor, uiReserveMinimum); }
    }
    void ResetMemoryReserveInfoToDefault()
    {
        if (m_priReserveInfo) { delete m_priReserveInfo; m_priReserveInfo = NULL; }
    }

    const dxWorldProcessMemoryManager *GetMemoryManager() const { return m_pmmMemoryManager; }
    const dxWorldProcessMemoryManager *SureGetMemoryManager() const { return m_pmmMemoryManager ? m_pmmMemoryManager : &g_WorldProcessMallocMemoryManager; }
    void SetMemoryManager(dxWorldProcessMemoryManager::alloc_block_fn_t fnAlloc, 
        dxWorldProcessMemoryManager::shrink_block_fn_t fnShrink, 
        dxWorldProcessMemoryManager::free_block_fn_t fnFree) 
    {
        if (m_pmmMemoryManager) { m_pmmMemoryManager->Assign(fnAlloc, fnShrink, fnFree); }
        else { m_pmmMemoryManager = new dxWorldProcessMemoryManager(fnAlloc, fnShrink, fnFree); }
    }
    void ResetMemoryManagerToDefault()
    {
        if (m_pmmMemoryManager) { delete m_pmmMemoryManager; m_pmmMemoryManager = NULL; }
    }

private:
    unsigned m_uiRefCount;
    dxWorldProcessContext *m_ppcProcessingContext;
    dxWorldProcessMemoryReserveInfo *m_priReserveInfo;
    dxWorldProcessMemoryManager *m_pmmMemoryManager;
};


#endif
