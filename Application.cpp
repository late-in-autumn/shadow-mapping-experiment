// Application.cpp: implementation of the Application class.
//
//////////////////////////////////////////////////////////////////////

#include "stdafx.h"
#include "CS580HW.h"
#include "Application.h"

#ifdef _DEBUG
#undef THIS_FILE
static char THIS_FILE[] = __FILE__;
#define new DEBUG_NEW
#endif

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

Application::Application()
{
	m_pRender = nullptr;		// the renderer
	m_pUserInput = nullptr;
	m_pFrameBuffer = nullptr;
}

Application::~Application()
{
	if (m_pUserInput != nullptr)
		delete m_pUserInput;
}
