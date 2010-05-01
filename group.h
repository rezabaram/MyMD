#ifndef GROUP_H
#define GROUP_H


///BROKEN!!!! stopped coding in the middle


template<typename T>
class CList{
	
	// ++++ member class for node
	template<typename U>
	class CNode{
		public:
		CNode(const U &v, CNode<U> *p=NULL, CNode<U> *n=NULL):_value(v), _next(n), _previous(p){};
		~CNode();

		U value()const {return _value;}; 
		U &value(){return _value;};
		CNode<U> *next()const {return _next;};
		CNode<U> *previous()const {return _previous;};

		private:
		
		U _value;
		CNode<U> * _next;
		CNode<U> * _previous;
	};
	// ---- member class for node

	public:
	class iterator{
		iterator(T *p):_p(p){};
		T operator *(){return (*_p);}
		T operator ++(){
			//_p=_next;
			return (*_p);
			}
		private:
		T * _p;
		};

	CList():_begin(NULL), _end(NULL), _size(0){};
	~CList(){};//FIXME

	void push_back(const T value){
		CNode<T> *p=new CNode<T>(value);
		if(_begin==NULL){
			assert(_size==0);
			_begin=p;
			}
		_end=p;
		_size++;
		return;
		}

		
	private:
	CNode<T> * _begin;
	CNode<T> * _end;
	long int _size;
};



#endif
