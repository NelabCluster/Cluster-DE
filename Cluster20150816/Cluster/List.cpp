#include "Common.h"
#include "List.h"

void AddToTail(ListNode** pHead, int value)
{
	ListNode* pNew = new ListNode();
	pNew->m_nValue = value;
	pNew->m_pNext = NULL;

	if(*pHead == NULL)
		*pHead = pNew;
	else
	{
		ListNode* pNode = *pHead;
		while(pNode->m_pNext != NULL)
			pNode = pNode->m_pNext;
		pNode->m_pNext = pNew;
	}
}

void RemoveNode(ListNode** pHead, int value)
{

}

void PrintNode(ListNode *pHead)
{
	if (pHead == NULL)
		return;

	ListNode *p;
	p = pHead;
	while(p!=NULL)
	{
		printf("%d\t",p->m_nValue);
		p = p->m_pNext;
	}
	printf("\n");
}