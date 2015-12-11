//链表节点
struct ListNode
{
	int m_nValue;
	ListNode* m_pNext;
};

//往链表尾部添加节点
void AddToTail(ListNode** pHead, int value);
//在链表中找到第一个含有某值的节点并删除改节点
void RemoveNode(ListNode** pHead, int value);

void PrintNode(ListNode *pHead);