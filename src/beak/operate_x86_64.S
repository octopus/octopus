	.file	"operate_ri_vec.c"
	.text
	.p2align 4,,15
.globl zoperate_as_
	.type	zoperate_as_, @function
zoperate_as_:
.LFB528:
	.cfi_startproc
	pushq	%r15
	.cfi_def_cfa_offset 16
	pushq	%r14
	.cfi_def_cfa_offset 24
	movq	%rdx, %r14
	.cfi_offset 14, -24
	.cfi_offset 15, -16
	pushq	%r13
	.cfi_def_cfa_offset 32
	pushq	%r12
	.cfi_def_cfa_offset 40
	movq	%r9, %r12
	.cfi_offset 12, -40
	.cfi_offset 13, -32
	pushq	%rbp
	.cfi_def_cfa_offset 48
	pushq	%rbx
	.cfi_def_cfa_offset 56
	subq	$9552, %rsp
	.cfi_def_cfa_offset 9608
	movl	(%rdi), %ebp
	.cfi_offset 3, -56
	.cfi_offset 6, -48
	movl	4(%rdi), %edi
	movq	%rcx, -72(%rsp)
	movq	9608(%rsp), %r13
	testl	%ebp, %ebp
	movl	%edi, -104(%rsp)
	jle	.L2
	leal	-1(%rbp), %ecx
	movq	(%rsi), %rbx
	leaq	-56(%rsp), %rdx
	xorl	%eax, %eax
	leaq	8(,%rcx,8), %rdi
	andl	$7, %ecx
	movq	%rbx, (%rdx,%rax,2)
	movq	%rbx, 8(%rdx,%rax,2)
	movb	$8, %al
	cmpq	%rdi, %rax
	je	.L2
	testq	%rcx, %rcx
	je	.L3
	cmpq	$1, %rcx
	je	.L171
	cmpq	$2, %rcx
	.p2align 4,,5
	.p2align 3
	je	.L172
	cmpq	$3, %rcx
	.p2align 4,,5
	.p2align 3
	je	.L173
	cmpq	$4, %rcx
	.p2align 4,,5
	.p2align 3
	je	.L174
	cmpq	$5, %rcx
	.p2align 4,,5
	.p2align 3
	je	.L175
	cmpq	$6, %rcx
	.p2align 4,,5
	.p2align 3
	je	.L176
	movq	8(%rsi), %r15
	movq	%r15, (%rdx,%rax,2)
	movq	%r15, 24(%rdx)
	movl	$16, %eax
.L176:
	movq	(%rsi,%rax), %rcx
	movq	%rcx, (%rdx,%rax,2)
	movq	%rcx, 8(%rdx,%rax,2)
	addq	$8, %rax
.L175:
	movq	(%rsi,%rax), %r11
	movq	%r11, (%rdx,%rax,2)
	movq	%r11, 8(%rdx,%rax,2)
	addq	$8, %rax
.L174:
	movq	(%rsi,%rax), %r10
	movq	%r10, (%rdx,%rax,2)
	movq	%r10, 8(%rdx,%rax,2)
	addq	$8, %rax
.L173:
	movq	(%rsi,%rax), %rbx
	movq	%rbx, (%rdx,%rax,2)
	movq	%rbx, 8(%rdx,%rax,2)
	addq	$8, %rax
.L172:
	movq	(%rsi,%rax), %r9
	movq	%r9, (%rdx,%rax,2)
	movq	%r9, 8(%rdx,%rax,2)
	addq	$8, %rax
.L171:
	movq	(%rsi,%rax), %r15
	movq	%r15, (%rdx,%rax,2)
	movq	%r15, 8(%rdx,%rax,2)
	addq	$8, %rax
	cmpq	%rdi, %rax
	je	.L2
.L3:
	movq	(%rsi,%rax), %r10
	leaq	8(%rax), %rcx
	leaq	16(%rax), %r15
	movq	%r10, (%rdx,%rax,2)
	movq	%r10, 8(%rdx,%rax,2)
	leaq	24(%rax), %r10
	movq	8(%rsi,%rax), %r9
	movq	%r9, (%rdx,%rcx,2)
	movq	%r9, 8(%rdx,%rcx,2)
	leaq	32(%rax), %rcx
	movq	16(%rsi,%rax), %rbx
	movq	%rbx, (%rdx,%r15,2)
	movq	%rbx, 8(%rdx,%r15,2)
	leaq	40(%rax), %r15
	movq	24(%rsi,%rax), %r11
	movq	%r11, (%rdx,%r10,2)
	movq	%r11, 8(%rdx,%r10,2)
	leaq	48(%rax), %r10
	movq	32(%rsi,%rax), %r9
	movq	%r9, (%rdx,%rcx,2)
	movq	%r9, 8(%rdx,%rcx,2)
	leaq	56(%rax), %rcx
	movq	40(%rsi,%rax), %rbx
	movq	%rbx, (%rdx,%r15,2)
	movq	%rbx, 8(%rdx,%r15,2)
	movq	48(%rsi,%rax), %r11
	movq	%r11, (%rdx,%r10,2)
	movq	%r11, 8(%rdx,%r10,2)
	movq	56(%rsi,%rax), %r9
	addq	$64, %rax
	cmpq	%rdi, %rax
	movq	%r9, (%rdx,%rcx,2)
	movq	%r9, 8(%rdx,%rcx,2)
	jne	.L3
.L2:
	testb	$15, %r12b
	jne	.L4
	testb	$15, %r13b
	je	.L5
.L4:
	cmpl	$0, -104(%rsp)
	.p2align 4,,4
	.p2align 3
	jle	.L36
	movslq	%ebp,%r11
	leaq	-56(%rsp), %rdx
	leaq	6344(%rsp), %rsi
	salq	$2, %r11
	movq	$0, -112(%rsp)
	movl	$0, -92(%rsp)
	movq	%r11, -88(%rsp)
.L35:
	movq	-72(%rsp), %rbx
	movq	-112(%rsp), %rax
	testl	%ebp, %ebp
	movl	(%rbx,%rax), %r15d
	jle	.L196
	movslq	(%r14),%rdi
	movslq	%r15d,%rcx
	leal	-1(%rbp), %ebx
	salq	$4, %rcx
	movl	$4, %eax
	andl	$3, %ebx
	salq	$4, %rdi
	leaq	(%r12,%rdi), %r9
	xorl	%edi, %edi
	cmpl	$1, %ebp
	movq	%r9, (%rsi,%rdi,2)
	movb	$1, %dil
	movupd	(%r9,%rcx), %xmm0
	mulpd	(%rdx), %xmm0
	jle	.L23
	testl	%ebx, %ebx
	je	.L24
	cmpl	$1, %ebx
	je	.L169
	cmpl	$2, %ebx
	.p2align 4,,5
	.p2align 3
	je	.L170
	movslq	4(%r14),%rax
	movl	$2, %edi
	salq	$4, %rax
	leaq	(%r12,%rax), %r9
	movl	$8, %eax
	movq	%r9, 8(%rsi)
	movupd	(%r9,%rcx), %xmm1
	mulpd	16(%rdx), %xmm1
	addpd	%xmm1, %xmm0
.L170:
	movslq	(%r14,%rax),%r10
	incl	%edi
	salq	$4, %r10
	leaq	(%r12,%r10), %r11
	movq	%r11, (%rsi,%rax,2)
	movupd	(%r11,%rcx), %xmm3
	mulpd	(%rdx,%rax,4), %xmm3
	addq	$4, %rax
	addpd	%xmm3, %xmm0
.L169:
	movslq	(%r14,%rax),%r9
	incl	%edi
	salq	$4, %r9
	leaq	(%r12,%r9), %rbx
	movq	%rbx, (%rsi,%rax,2)
	movupd	(%rbx,%rcx), %xmm2
	mulpd	(%rdx,%rax,4), %xmm2
	addq	$4, %rax
	cmpl	%edi, %ebp
	addpd	%xmm2, %xmm0
	jle	.L23
	.p2align 4,,7
	.p2align 3
.L24:
	movslq	(%r14,%rax),%r9
	movslq	4(%r14,%rax),%r10
	addl	$4, %edi
	salq	$4, %r9
	salq	$4, %r10
	leaq	(%r12,%r9), %r11
	leaq	(%r12,%r10), %rbx
	leaq	4(%rax), %r9
	movq	%r11, (%rsi,%rax,2)
	movupd	(%r11,%rcx), %xmm6
	movslq	8(%r14,%rax),%r11
	movq	%rbx, (%rsi,%r9,2)
	mulpd	(%rdx,%rax,4), %xmm6
	movupd	(%rbx,%rcx), %xmm7
	movslq	12(%r14,%rax),%rbx
	salq	$4, %r11
	mulpd	(%rdx,%r9,4), %xmm7
	leaq	(%r12,%r11), %r10
	leaq	8(%rax), %r9
	movq	%r10, (%rsi,%r9,2)
	salq	$4, %rbx
	addpd	%xmm6, %xmm0
	movupd	(%r10,%rcx), %xmm2
	leaq	(%r12,%rbx), %r11
	leaq	12(%rax), %r10
	addq	$16, %rax
	mulpd	(%rdx,%r9,4), %xmm2
	cmpl	%edi, %ebp
	movq	%r11, (%rsi,%r10,2)
	movupd	(%r11,%rcx), %xmm4
	addpd	%xmm7, %xmm0
	mulpd	(%rdx,%r10,4), %xmm4
	addpd	%xmm2, %xmm0
	addpd	%xmm4, %xmm0
	jg	.L24
.L23:
	movupd	%xmm0, (%r13,%rcx)
	movq	-112(%rsp), %rdi
	incl	%r15d
	movl	(%r8,%rdi), %eax
	leal	-3(%rax), %ecx
	cmpl	%ecx, %r15d
	jge	.L25
	movslq	%r15d,%rdi
	movq	%r14, -80(%rsp)
	movq	-112(%rsp), %r14
	leaq	1(%rdi), %r10
	leaq	2(%rdi), %r9
	movq	%rdi, %r11
	addq	$3, %rdi
	salq	$4, %r11
	salq	$4, %r10
	salq	$4, %r9
	salq	$4, %rdi
	.p2align 4,,7
	.p2align 3
.L26:
	testl	%ebp, %ebp
	jle	.L197
	movq	(%rsi), %rax
	movlpd	(%rdx), %xmm8
	leal	-1(%rbp), %ecx
	movl	$1, %ebx
	movhpd	8(%rdx), %xmm8
	andl	$1, %ecx
	cmpl	$1, %ebp
	movupd	(%rax,%r11), %xmm2
	movupd	(%rax,%r10), %xmm3
	mulpd	%xmm8, %xmm2
	movupd	(%rax,%r9), %xmm4
	mulpd	%xmm8, %xmm3
	movupd	(%rax,%rdi), %xmm1
	movl	$8, %eax
	mulpd	%xmm8, %xmm4
	mulpd	%xmm8, %xmm1
	jle	.L29
	testl	%ecx, %ecx
	je	.L27
	movq	8(%rsi), %rbx
	movlpd	16(%rdx), %xmm14
	cmpl	$2, %ebp
	movl	$16, %eax
	movhpd	24(%rdx), %xmm14
	movupd	(%rbx,%r11), %xmm13
	movupd	(%rbx,%r10), %xmm12
	mulpd	%xmm14, %xmm13
	movupd	(%rbx,%r9), %xmm11
	mulpd	%xmm14, %xmm12
	movupd	(%rbx,%rdi), %xmm10
	movl	$2, %ebx
	mulpd	%xmm14, %xmm11
	addpd	%xmm13, %xmm2
	mulpd	%xmm14, %xmm10
	addpd	%xmm12, %xmm3
	addpd	%xmm11, %xmm4
	addpd	%xmm10, %xmm1
	jle	.L29
	.p2align 4,,7
	.p2align 3
.L27:
	movq	(%rsi,%rax), %rcx
	movlpd	(%rdx,%rax,2), %xmm8
	addl	$2, %ebx
	movlpd	16(%rdx,%rax,2), %xmm13
	movhpd	8(%rdx,%rax,2), %xmm8
	movupd	(%rcx,%r11), %xmm0
	movhpd	24(%rdx,%rax,2), %xmm13
	movupd	(%rcx,%r10), %xmm6
	mulpd	%xmm8, %xmm0
	movupd	(%rcx,%r9), %xmm7
	mulpd	%xmm8, %xmm6
	movupd	(%rcx,%rdi), %xmm14
	mulpd	%xmm8, %xmm7
	movq	8(%rsi,%rax), %rcx
	addq	$16, %rax
	addpd	%xmm0, %xmm2
	cmpl	%ebx, %ebp
	mulpd	%xmm8, %xmm14
	movupd	(%rcx,%r11), %xmm12
	addpd	%xmm6, %xmm3
	movupd	(%rcx,%r10), %xmm11
	addpd	%xmm7, %xmm4
	mulpd	%xmm13, %xmm12
	movupd	(%rcx,%r9), %xmm10
	addpd	%xmm14, %xmm1
	mulpd	%xmm13, %xmm11
	movupd	(%rcx,%rdi), %xmm5
	mulpd	%xmm13, %xmm10
	addpd	%xmm12, %xmm2
	mulpd	%xmm13, %xmm5
	addpd	%xmm11, %xmm3
	addpd	%xmm10, %xmm4
	addpd	%xmm5, %xmm1
	jg	.L27
.L29:
	movupd	%xmm2, (%r13,%r11)
	addl	$4, %r15d
	addq	$64, %r11
	movupd	%xmm3, (%r13,%r10)
	addq	$64, %r10
	movupd	%xmm4, (%r13,%r9)
	addq	$64, %r9
	movupd	%xmm1, (%r13,%rdi)
	addq	$64, %rdi
	movl	(%r8,%r14), %eax
	leal	-3(%rax), %ebx
	cmpl	%r15d, %ebx
	jg	.L26
	movq	-80(%rsp), %r14
.L25:
	cmpl	%eax, %r15d
	jge	.L30
	movq	-112(%rsp), %r9
	movslq	%r15d,%rcx
	salq	$4, %rcx
	.p2align 4,,7
	.p2align 3
.L31:
	testl	%ebp, %ebp
	jle	.L198
	movq	(%rsi), %r10
	leal	-1(%rbp), %edi
	movl	$1, %ebx
	movl	$8, %eax
	andl	$7, %edi
	addq	%rcx, %r10
	cmpl	$1, %ebp
	movupd	(%r10), %xmm0
	mulpd	(%rdx), %xmm0
	jle	.L34
	testl	%edi, %edi
	je	.L32
	cmpl	$1, %edi
	je	.L163
	cmpl	$2, %edi
	.p2align 4,,5
	.p2align 3
	je	.L164
	cmpl	$3, %edi
	.p2align 4,,5
	.p2align 3
	je	.L165
	cmpl	$4, %edi
	.p2align 4,,5
	.p2align 3
	je	.L166
	cmpl	$5, %edi
	.p2align 4,,5
	.p2align 3
	je	.L167
	cmpl	$6, %edi
	.p2align 4,,5
	.p2align 3
	je	.L168
	movq	8(%rsi), %rax
	movl	$2, %ebx
	addq	%rcx, %rax
	movupd	(%rax), %xmm15
	movl	$16, %eax
	mulpd	16(%rdx), %xmm15
	addpd	%xmm15, %xmm0
.L168:
	movq	(%rsi,%rax), %r11
	incl	%ebx
	addq	%rcx, %r11
	movupd	(%r11), %xmm7
	mulpd	(%rdx,%rax,2), %xmm7
	addq	$8, %rax
	addpd	%xmm7, %xmm0
.L167:
	movq	(%rsi,%rax), %rdi
	incl	%ebx
	addq	%rcx, %rdi
	movupd	(%rdi), %xmm6
	mulpd	(%rdx,%rax,2), %xmm6
	addq	$8, %rax
	addpd	%xmm6, %xmm0
.L166:
	movq	(%rsi,%rax), %r10
	incl	%ebx
	addq	%rcx, %r10
	movupd	(%r10), %xmm8
	mulpd	(%rdx,%rax,2), %xmm8
	addq	$8, %rax
	addpd	%xmm8, %xmm0
.L165:
	movq	(%rsi,%rax), %r11
	incl	%ebx
	addq	%rcx, %r11
	movupd	(%r11), %xmm4
	mulpd	(%rdx,%rax,2), %xmm4
	addq	$8, %rax
	addpd	%xmm4, %xmm0
.L164:
	movq	(%rsi,%rax), %rdi
	incl	%ebx
	addq	%rcx, %rdi
	movupd	(%rdi), %xmm9
	mulpd	(%rdx,%rax,2), %xmm9
	addq	$8, %rax
	addpd	%xmm9, %xmm0
.L163:
	movq	(%rsi,%rax), %r10
	incl	%ebx
	addq	%rcx, %r10
	movupd	(%r10), %xmm5
	mulpd	(%rdx,%rax,2), %xmm5
	addq	$8, %rax
	cmpl	%ebx, %ebp
	addpd	%xmm5, %xmm0
	jle	.L34
	.p2align 4,,7
	.p2align 3
.L32:
	movq	(%rsi,%rax), %rdi
	movq	8(%rsi,%rax), %r11
	addl	$8, %ebx
	movq	16(%rsi,%rax), %r10
	addq	%rcx, %rdi
	addq	%rcx, %r11
	movupd	(%rdi), %xmm11
	addq	%rcx, %r10
	movq	24(%rsi,%rax), %rdi
	mulpd	(%rdx,%rax,2), %xmm11
	movupd	(%r11), %xmm10
	addq	%rcx, %rdi
	movq	32(%rsi,%rax), %r11
	mulpd	16(%rdx,%rax,2), %xmm10
	movupd	(%r10), %xmm9
	addq	%rcx, %r11
	movq	40(%rsi,%rax), %r10
	mulpd	32(%rdx,%rax,2), %xmm9
	addpd	%xmm11, %xmm0
	movupd	(%rdi), %xmm5
	addq	%rcx, %r10
	movq	48(%rsi,%rax), %rdi
	mulpd	48(%rdx,%rax,2), %xmm5
	movupd	(%r11), %xmm2
	addq	%rcx, %rdi
	movq	56(%rsi,%rax), %r11
	addpd	%xmm10, %xmm0
	mulpd	64(%rdx,%rax,2), %xmm2
	movupd	(%r10), %xmm4
	addq	%rcx, %r11
	mulpd	80(%rdx,%rax,2), %xmm4
	movupd	(%rdi), %xmm3
	addpd	%xmm9, %xmm0
	mulpd	96(%rdx,%rax,2), %xmm3
	movupd	(%r11), %xmm1
	mulpd	112(%rdx,%rax,2), %xmm1
	addq	$64, %rax
	addpd	%xmm5, %xmm0
	cmpl	%ebx, %ebp
	addpd	%xmm2, %xmm0
	addpd	%xmm4, %xmm0
	addpd	%xmm3, %xmm0
	addpd	%xmm1, %xmm0
	jg	.L32
.L34:
	movupd	%xmm0, (%r13,%rcx)
	incl	%r15d
	addq	$16, %rcx
	cmpl	%r15d, (%r8,%r9)
	jg	.L31
.L30:
	incl	-92(%rsp)
	addq	$4, -112(%rsp)
	movl	-92(%rsp), %r15d
	addq	-88(%rsp), %r14
	cmpl	%r15d, -104(%rsp)
	jg	.L35
.L36:
	addq	$9552, %rsp
	popq	%rbx
	popq	%rbp
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	ret
.L5:
	cmpl	$0, -104(%rsp)
	jle	.L36
	movl	-104(%rsp), %esi
	leal	-1(%rbp), %eax
	movslq	%ebp,%rdi
	salq	$2, %rdi
	movq	$0, -104(%rsp)
	incq	%rax
	movq	%rdi, -64(%rsp)
	movq	%rax, -88(%rsp)
	decl	%esi
	movq	-88(%rsp), %r11
	leaq	4(,%rsi,4), %rdx
	leaq	-56(%rsp), %rax
	movq	%rdx, -80(%rsp)
	salq	$3, %r11
.L21:
	movq	-72(%rsp), %rcx
	movq	-104(%rsp), %rbx
	testl	%ebp, %ebp
	movl	(%rcx,%rbx), %r15d
	jle	.L199
	movslq	(%r14),%rdx
	movq	-88(%rsp), %r9
	movslq	%r15d,%rcx
	movlpd	(%rax), %xmm1
	salq	$4, %rcx
	leaq	6344(%rsp), %rbx
	salq	$2, %r9
	movhpd	8(%rax), %xmm1
	salq	$4, %rdx
	leaq	-4(%r9), %rsi
	leaq	(%r12,%rdx), %r10
	movl	$4, %edx
	shrq	$2, %rsi
	movq	%r10, 6344(%rsp)
	andl	$3, %esi
	cmpq	%r9, %rdx
	movapd	(%r10,%rcx), %xmm0
	mulpd	%xmm1, %xmm0
	je	.L9
	testq	%rsi, %rsi
	je	.L10
	cmpq	$1, %rsi
	je	.L161
	cmpq	$2, %rsi
	.p2align 4,,5
	.p2align 3
	je	.L162
	movslq	4(%r14),%r10
	movlpd	16(%rax), %xmm3
	movl	$8, %edx
	movhpd	24(%rax), %xmm3
	salq	$4, %r10
	leaq	(%r12,%r10), %rdi
	movq	%rdi, 8(%rbx)
	movapd	(%rdi,%rcx), %xmm1
	mulpd	%xmm3, %xmm1
	addpd	%xmm1, %xmm0
.L162:
	movslq	(%r14,%rdx),%rdi
	movlpd	(%rax,%rdx,4), %xmm2
	movhpd	8(%rax,%rdx,4), %xmm2
	salq	$4, %rdi
	leaq	(%r12,%rdi), %rsi
	movq	%rsi, (%rbx,%rdx,2)
	addq	$4, %rdx
	movapd	(%rsi,%rcx), %xmm12
	mulpd	%xmm2, %xmm12
	addpd	%xmm12, %xmm0
.L161:
	movslq	(%r14,%rdx),%rsi
	movlpd	(%rax,%rdx,4), %xmm14
	movhpd	8(%rax,%rdx,4), %xmm14
	salq	$4, %rsi
	leaq	(%r12,%rsi), %r10
	movq	%r10, (%rbx,%rdx,2)
	addq	$4, %rdx
	movapd	(%r10,%rcx), %xmm13
	cmpq	%r9, %rdx
	mulpd	%xmm14, %xmm13
	addpd	%xmm13, %xmm0
	je	.L9
	.p2align 4,,7
	.p2align 3
.L10:
	movslq	(%r14,%rdx),%rsi
	movslq	4(%r14,%rdx),%r10
	movlpd	(%rax,%rdx,4), %xmm9
	movhpd	8(%rax,%rdx,4), %xmm9
	salq	$4, %rsi
	salq	$4, %r10
	leaq	(%r12,%rsi), %rdi
	leaq	4(%rdx), %rsi
	movq	%rdi, (%rbx,%rdx,2)
	movlpd	(%rax,%rsi,4), %xmm5
	movapd	(%rdi,%rcx), %xmm8
	leaq	(%r12,%r10), %rdi
	movslq	8(%r14,%rdx),%r10
	movhpd	8(%rax,%rsi,4), %xmm5
	mulpd	%xmm9, %xmm8
	movq	%rdi, (%rbx,%rsi,2)
	leaq	8(%rdx), %rsi
	movapd	(%rdi,%rcx), %xmm7
	salq	$4, %r10
	mulpd	%xmm5, %xmm7
	leaq	(%r12,%r10), %rdi
	movslq	12(%r14,%rdx),%r10
	movlpd	(%rax,%rsi,4), %xmm6
	addpd	%xmm8, %xmm0
	movq	%rdi, (%rbx,%rsi,2)
	movhpd	8(%rax,%rsi,4), %xmm6
	leaq	12(%rdx), %rsi
	movapd	(%rdi,%rcx), %xmm4
	addq	$16, %rdx
	salq	$4, %r10
	cmpq	%r9, %rdx
	leaq	(%r12,%r10), %rdi
	mulpd	%xmm6, %xmm4
	addpd	%xmm7, %xmm0
	movlpd	(%rax,%rsi,4), %xmm2
	movq	%rdi, (%rbx,%rsi,2)
	movhpd	8(%rax,%rsi,4), %xmm2
	movapd	(%rdi,%rcx), %xmm3
	mulpd	%xmm2, %xmm3
	addpd	%xmm4, %xmm0
	addpd	%xmm3, %xmm0
	jne	.L10
.L9:
	movq	-104(%rsp), %r9
	incl	%r15d
	movaps	%xmm0, (%r13,%rcx)
	movl	(%r8,%r9), %edx
	movl	%edx, -92(%rsp)
	subl	$3, %edx
	cmpl	%edx, %r15d
	jge	.L11
	movslq	%r15d,%rbx
	leaq	6344(%rsp), %r10
	movl	%edx, -112(%rsp)
	leaq	1(%rbx), %rdi
	leaq	2(%rbx), %rsi
	movq	%rbx, %r9
	addq	$3, %rbx
	salq	$4, %r9
	salq	$4, %rdi
	salq	$4, %rsi
	salq	$4, %rbx
	.p2align 4,,7
	.p2align 3
.L12:
	testl	%ebp, %ebp
	jle	.L200
	movlpd	(%rax), %xmm1
	movq	(%r10), %rdx
	leaq	-8(%r11), %rcx
	movhpd	8(%rax), %xmm1
	shrq	$3, %rcx
	movapd	(%rdx,%rdi), %xmm3
	andl	$3, %ecx
	movapd	%xmm1, %xmm0
	mulpd	(%rdx,%r9), %xmm1
	movapd	(%rdx,%rsi), %xmm2
	mulpd	%xmm0, %xmm3
	mulpd	%xmm0, %xmm2
	mulpd	(%rdx,%rbx), %xmm0
	movl	$8, %edx
	cmpq	%r11, %rdx
	je	.L15
	testq	%rcx, %rcx
	je	.L13
	cmpq	$1, %rcx
	je	.L159
	cmpq	$2, %rcx
	.p2align 4,,5
	.p2align 3
	je	.L160
	movlpd	16(%rax), %xmm14
	movq	8(%r10), %rdx
	movhpd	24(%rax), %xmm14
	movapd	(%rdx,%rdi), %xmm13
	movapd	%xmm14, %xmm11
	mulpd	(%rdx,%r9), %xmm14
	movapd	(%rdx,%rsi), %xmm12
	mulpd	%xmm11, %xmm13
	mulpd	%xmm11, %xmm12
	mulpd	(%rdx,%rbx), %xmm11
	movl	$16, %edx
	addpd	%xmm14, %xmm1
	addpd	%xmm13, %xmm3
	addpd	%xmm12, %xmm2
	addpd	%xmm11, %xmm0
.L160:
	movlpd	(%rax,%rdx,2), %xmm4
	movq	(%r10,%rdx), %rcx
	movhpd	8(%rax,%rdx,2), %xmm4
	addq	$8, %rdx
	movapd	(%rcx,%r9), %xmm8
	movapd	%xmm4, %xmm15
	movapd	(%rcx,%rdi), %xmm6
	mulpd	%xmm4, %xmm8
	movapd	(%rcx,%rsi), %xmm7
	mulpd	(%rcx,%rbx), %xmm15
	mulpd	%xmm4, %xmm6
	addpd	%xmm8, %xmm1
	mulpd	%xmm4, %xmm7
	addpd	%xmm6, %xmm3
	addpd	%xmm7, %xmm2
	addpd	%xmm15, %xmm0
.L159:
	movlpd	(%rax,%rdx,2), %xmm11
	movq	(%r10,%rdx), %rcx
	movhpd	8(%rax,%rdx,2), %xmm11
	addq	$8, %rdx
	movapd	(%rcx,%rdi), %xmm10
	cmpq	%r11, %rdx
	movapd	%xmm11, %xmm9
	mulpd	(%rcx,%r9), %xmm11
	movapd	(%rcx,%rsi), %xmm5
	mulpd	%xmm9, %xmm10
	mulpd	%xmm9, %xmm5
	mulpd	(%rcx,%rbx), %xmm9
	addpd	%xmm11, %xmm1
	addpd	%xmm10, %xmm3
	addpd	%xmm5, %xmm2
	addpd	%xmm9, %xmm0
	je	.L15
	.p2align 4,,7
	.p2align 3
.L13:
	movlpd	(%rax,%rdx,2), %xmm15
	movq	(%r10,%rdx), %rcx
	movlpd	16(%rax,%rdx,2), %xmm11
	movhpd	8(%rax,%rdx,2), %xmm15
	movlpd	32(%rax,%rdx,2), %xmm7
	movapd	(%rcx,%rdi), %xmm14
	movhpd	24(%rax,%rdx,2), %xmm11
	movapd	%xmm15, %xmm12
	mulpd	(%rcx,%r9), %xmm15
	movapd	(%rcx,%rsi), %xmm13
	movhpd	40(%rax,%rdx,2), %xmm7
	mulpd	%xmm12, %xmm14
	movapd	%xmm11, %xmm8
	mulpd	%xmm12, %xmm13
	mulpd	(%rcx,%rbx), %xmm12
	movq	8(%r10,%rdx), %rcx
	addpd	%xmm15, %xmm1
	movapd	%xmm7, %xmm15
	addpd	%xmm14, %xmm3
	movlpd	48(%rax,%rdx,2), %xmm14
	movapd	(%rcx,%r9), %xmm10
	mulpd	(%rcx,%rbx), %xmm8
	addpd	%xmm13, %xmm2
	movapd	(%rcx,%rdi), %xmm9
	mulpd	%xmm11, %xmm10
	addpd	%xmm12, %xmm0
	movapd	(%rcx,%rsi), %xmm5
	movq	16(%r10,%rdx), %rcx
	mulpd	%xmm11, %xmm9
	movhpd	56(%rax,%rdx,2), %xmm14
	mulpd	%xmm11, %xmm5
	movapd	(%rcx,%rdi), %xmm6
	addpd	%xmm10, %xmm1
	mulpd	(%rcx,%r9), %xmm7
	movapd	(%rcx,%rsi), %xmm4
	addpd	%xmm9, %xmm3
	mulpd	%xmm15, %xmm6
	addpd	%xmm5, %xmm2
	mulpd	%xmm15, %xmm4
	addpd	%xmm8, %xmm0
	mulpd	(%rcx,%rbx), %xmm15
	movq	24(%r10,%rdx), %rcx
	movapd	%xmm14, %xmm10
	addpd	%xmm7, %xmm1
	addq	$32, %rdx
	cmpq	%r11, %rdx
	movapd	(%rcx,%r9), %xmm13
	addpd	%xmm6, %xmm3
	mulpd	(%rcx,%rbx), %xmm10
	movapd	(%rcx,%rdi), %xmm12
	addpd	%xmm4, %xmm2
	mulpd	%xmm14, %xmm13
	movapd	(%rcx,%rsi), %xmm11
	addpd	%xmm15, %xmm0
	mulpd	%xmm14, %xmm12
	mulpd	%xmm14, %xmm11
	addpd	%xmm13, %xmm1
	addpd	%xmm12, %xmm3
	addpd	%xmm11, %xmm2
	addpd	%xmm10, %xmm0
	jne	.L13
.L15:
	movaps	%xmm1, (%r13,%r9)
	addl	$4, %r15d
	addq	$64, %r9
	movaps	%xmm3, (%r13,%rdi)
	addq	$64, %rdi
	movaps	%xmm2, (%r13,%rsi)
	addq	$64, %rsi
	movaps	%xmm0, (%r13,%rbx)
	addq	$64, %rbx
	cmpl	%r15d, -112(%rsp)
	jg	.L12
.L11:
	cmpl	-92(%rsp), %r15d
	jge	.L16
	movq	-104(%rsp), %rdi
	movslq	%r15d,%rcx
	leaq	6344(%rsp), %rbx
	salq	$4, %rcx
	.p2align 4,,7
	.p2align 3
.L17:
	testl	%ebp, %ebp
	jle	.L201
	movq	(%rbx), %rdx
	leaq	-8(%r11), %rsi
	movlpd	(%rax), %xmm3
	shrq	$3, %rsi
	movhpd	8(%rax), %xmm3
	movapd	(%rdx,%rcx), %xmm0
	movl	$8, %edx
	andl	$7, %esi
	cmpq	%r11, %rdx
	mulpd	%xmm3, %xmm0
	je	.L20
	testq	%rsi, %rsi
	je	.L18
	cmpq	$1, %rsi
	je	.L153
	cmpq	$2, %rsi
	.p2align 4,,5
	.p2align 3
	je	.L154
	cmpq	$3, %rsi
	.p2align 4,,5
	.p2align 3
	je	.L155
	cmpq	$4, %rsi
	.p2align 4,,5
	.p2align 3
	je	.L156
	cmpq	$5, %rsi
	.p2align 4,,5
	.p2align 3
	je	.L157
	cmpq	$6, %rsi
	.p2align 4,,5
	.p2align 3
	je	.L158
	movq	8(%rbx), %rsi
	movlpd	16(%rax), %xmm13
	movl	$16, %edx
	movhpd	24(%rax), %xmm13
	movapd	(%rsi,%rcx), %xmm12
	mulpd	%xmm13, %xmm12
	addpd	%xmm12, %xmm0
.L158:
	movq	(%rbx,%rdx), %r9
	movlpd	(%rax,%rdx,2), %xmm15
	movhpd	8(%rax,%rdx,2), %xmm15
	addq	$8, %rdx
	movapd	(%r9,%rcx), %xmm14
	mulpd	%xmm15, %xmm14
	addpd	%xmm14, %xmm0
.L157:
	movq	(%rbx,%rdx), %r10
	movlpd	(%rax,%rdx,2), %xmm6
	movhpd	8(%rax,%rdx,2), %xmm6
	addq	$8, %rdx
	movapd	(%r10,%rcx), %xmm7
	mulpd	%xmm6, %xmm7
	addpd	%xmm7, %xmm0
.L156:
	movq	(%rbx,%rdx), %rsi
	movlpd	(%rax,%rdx,2), %xmm1
	movhpd	8(%rax,%rdx,2), %xmm1
	addq	$8, %rdx
	movapd	(%rsi,%rcx), %xmm8
	mulpd	%xmm1, %xmm8
	addpd	%xmm8, %xmm0
.L155:
	movq	(%rbx,%rdx), %r9
	movlpd	(%rax,%rdx,2), %xmm4
	movhpd	8(%rax,%rdx,2), %xmm4
	addq	$8, %rdx
	movapd	(%r9,%rcx), %xmm3
	mulpd	%xmm4, %xmm3
	addpd	%xmm3, %xmm0
.L154:
	movq	(%rbx,%rdx), %r10
	movlpd	(%rax,%rdx,2), %xmm5
	movhpd	8(%rax,%rdx,2), %xmm5
	addq	$8, %rdx
	movapd	(%r10,%rcx), %xmm2
	mulpd	%xmm5, %xmm2
	addpd	%xmm2, %xmm0
.L153:
	movq	(%rbx,%rdx), %rsi
	movlpd	(%rax,%rdx,2), %xmm10
	movhpd	8(%rax,%rdx,2), %xmm10
	addq	$8, %rdx
	movapd	(%rsi,%rcx), %xmm9
	cmpq	%r11, %rdx
	mulpd	%xmm10, %xmm9
	addpd	%xmm9, %xmm0
	je	.L20
	.p2align 4,,7
	.p2align 3
.L18:
	movq	(%rbx,%rdx), %rsi
	movlpd	(%rax,%rdx,2), %xmm1
	movq	8(%rbx,%rdx), %r10
	movlpd	16(%rax,%rdx,2), %xmm14
	movhpd	8(%rax,%rdx,2), %xmm1
	movq	16(%rbx,%rdx), %r9
	movlpd	32(%rax,%rdx,2), %xmm12
	movapd	(%rsi,%rcx), %xmm15
	movhpd	24(%rax,%rdx,2), %xmm14
	movhpd	40(%rax,%rdx,2), %xmm12
	movq	24(%rbx,%rdx), %rsi
	mulpd	%xmm1, %xmm15
	movapd	(%r10,%rcx), %xmm13
	movq	32(%rbx,%rdx), %r10
	mulpd	%xmm14, %xmm13
	movapd	(%r9,%rcx), %xmm11
	movq	40(%rbx,%rdx), %r9
	mulpd	%xmm12, %xmm11
	movlpd	48(%rax,%rdx,2), %xmm10
	addpd	%xmm15, %xmm0
	movapd	(%rsi,%rcx), %xmm9
	movhpd	56(%rax,%rdx,2), %xmm10
	movq	48(%rbx,%rdx), %rsi
	movlpd	64(%rax,%rdx,2), %xmm5
	mulpd	%xmm10, %xmm9
	movapd	(%r10,%rcx), %xmm8
	movhpd	72(%rax,%rdx,2), %xmm5
	movq	56(%rbx,%rdx), %r10
	addpd	%xmm13, %xmm0
	movlpd	80(%rax,%rdx,2), %xmm6
	mulpd	%xmm5, %xmm8
	movapd	(%r9,%rcx), %xmm7
	movhpd	88(%rax,%rdx,2), %xmm6
	movlpd	96(%rax,%rdx,2), %xmm2
	mulpd	%xmm6, %xmm7
	movapd	(%rsi,%rcx), %xmm4
	addpd	%xmm11, %xmm0
	movhpd	104(%rax,%rdx,2), %xmm2
	movlpd	112(%rax,%rdx,2), %xmm3
	movapd	(%r10,%rcx), %xmm1
	mulpd	%xmm2, %xmm4
	movhpd	120(%rax,%rdx,2), %xmm3
	addq	$64, %rdx
	cmpq	%r11, %rdx
	addpd	%xmm9, %xmm0
	mulpd	%xmm3, %xmm1
	addpd	%xmm8, %xmm0
	addpd	%xmm7, %xmm0
	addpd	%xmm4, %xmm0
	addpd	%xmm1, %xmm0
	jne	.L18
.L20:
	movaps	%xmm0, (%r13,%rcx)
	incl	%r15d
	addq	$16, %rcx
	cmpl	%r15d, (%r8,%rdi)
	jg	.L17
.L16:
	addq	$4, -104(%rsp)
	addq	-64(%rsp), %r14
	movq	-80(%rsp), %r15
	cmpq	%r15, -104(%rsp)
	jne	.L21
	jmp	.L36
.L201:
	xorpd	%xmm0, %xmm0
	jmp	.L20
.L200:
	xorpd	%xmm0, %xmm0
	movapd	%xmm0, %xmm2
	movapd	%xmm0, %xmm3
	movapd	%xmm0, %xmm1
	jmp	.L15
.L198:
	xorpd	%xmm0, %xmm0
	jmp	.L34
.L197:
	xorpd	%xmm1, %xmm1
	movapd	%xmm1, %xmm4
	movapd	%xmm1, %xmm3
	movapd	%xmm1, %xmm2
	jmp	.L29
.L196:
	movslq	%r15d,%rcx
	xorpd	%xmm0, %xmm0
	salq	$4, %rcx
	jmp	.L23
.L199:
	movslq	%r15d,%rcx
	xorpd	%xmm0, %xmm0
	salq	$4, %rcx
	jmp	.L9
	.cfi_endproc
.LFE528:
	.size	zoperate_as_, .-zoperate_as_
	.p2align 4,,15
.globl doperate_as_
	.type	doperate_as_, @function
doperate_as_:
.LFB527:
	.cfi_startproc
	pushq	%r15
	.cfi_def_cfa_offset 16
	pushq	%r14
	.cfi_def_cfa_offset 24
	pushq	%r13
	.cfi_def_cfa_offset 32
	pushq	%r12
	.cfi_def_cfa_offset 40
	pushq	%rbp
	.cfi_def_cfa_offset 48
	pushq	%rbx
	.cfi_def_cfa_offset 56
	subq	$6352, %rsp
	.cfi_def_cfa_offset 6408
	movl	(%rdi), %r11d
	movl	4(%rdi), %edi
	movq	%rcx, -72(%rsp)
	movq	%r8, -88(%rsp)
	testl	%r11d, %r11d
	movl	%edi, -92(%rsp)
	jle	.L203
	.cfi_offset 3, -56
	.cfi_offset 6, -48
	.cfi_offset 12, -40
	.cfi_offset 13, -32
	.cfi_offset 14, -24
	.cfi_offset 15, -16
	leal	-1(%r11), %ebx
	movq	(%rsi), %rbp
	leaq	-56(%rsp), %rcx
	xorl	%eax, %eax
	leaq	8(,%rbx,8), %rdi
	andl	$7, %ebx
	movq	%rbp, (%rcx,%rax,2)
	movq	%rbp, 8(%rcx,%rax,2)
	movb	$8, %al
	cmpq	%rdi, %rax
	je	.L203
	testq	%rbx, %rbx
	je	.L204
	cmpq	$1, %rbx
	je	.L279
	cmpq	$2, %rbx
	.p2align 4,,5
	.p2align 3
	je	.L280
	cmpq	$3, %rbx
	.p2align 4,,5
	.p2align 3
	je	.L281
	cmpq	$4, %rbx
	.p2align 4,,5
	.p2align 3
	je	.L282
	cmpq	$5, %rbx
	.p2align 4,,5
	.p2align 3
	je	.L283
	cmpq	$6, %rbx
	.p2align 4,,5
	.p2align 3
	je	.L284
	movq	8(%rsi), %rbx
	movq	%rbx, (%rcx,%rax,2)
	movq	%rbx, 24(%rcx)
	movl	$16, %eax
.L284:
	movq	(%rsi,%rax), %r8
	movq	%r8, (%rcx,%rax,2)
	movq	%r8, 8(%rcx,%rax,2)
	addq	$8, %rax
.L283:
	movq	(%rsi,%rax), %r12
	movq	%r12, (%rcx,%rax,2)
	movq	%r12, 8(%rcx,%rax,2)
	addq	$8, %rax
.L282:
	movq	(%rsi,%rax), %r13
	movq	%r13, (%rcx,%rax,2)
	movq	%r13, 8(%rcx,%rax,2)
	addq	$8, %rax
.L281:
	movq	(%rsi,%rax), %r15
	movq	%r15, (%rcx,%rax,2)
	movq	%r15, 8(%rcx,%rax,2)
	addq	$8, %rax
.L280:
	movq	(%rsi,%rax), %r10
	movq	%r10, (%rcx,%rax,2)
	movq	%r10, 8(%rcx,%rax,2)
	addq	$8, %rax
.L279:
	movq	(%rsi,%rax), %rbp
	movq	%rbp, (%rcx,%rax,2)
	movq	%rbp, 8(%rcx,%rax,2)
	addq	$8, %rax
	cmpq	%rdi, %rax
	je	.L203
.L204:
	movq	(%rsi,%rax), %r15
	leaq	8(%rax), %r13
	leaq	16(%rax), %r10
	leaq	24(%rax), %rbx
	movq	%r15, (%rcx,%rax,2)
	movq	%r15, 8(%rcx,%rax,2)
	leaq	32(%rax), %r15
	movq	8(%rsi,%rax), %r14
	movq	%r14, (%rcx,%r13,2)
	movq	%r14, 8(%rcx,%r13,2)
	leaq	40(%rax), %r13
	movq	16(%rsi,%rax), %r12
	movq	%r12, (%rcx,%r10,2)
	movq	%r12, 8(%rcx,%r10,2)
	leaq	48(%rax), %r10
	movq	24(%rsi,%rax), %r8
	movq	%r8, (%rcx,%rbx,2)
	movq	%r8, 8(%rcx,%rbx,2)
	leaq	56(%rax), %rbx
	movq	32(%rsi,%rax), %rbp
	movq	%rbp, (%rcx,%r15,2)
	movq	%rbp, 8(%rcx,%r15,2)
	movq	40(%rsi,%rax), %r14
	movq	%r14, (%rcx,%r13,2)
	movq	%r14, 8(%rcx,%r13,2)
	movq	48(%rsi,%rax), %r12
	movq	%r12, (%rcx,%r10,2)
	movq	%r12, 8(%rcx,%r10,2)
	movq	56(%rsi,%rax), %r8
	addq	$64, %rax
	cmpq	%rdi, %rax
	movq	%r8, (%rcx,%rbx,2)
	movq	%r8, 8(%rcx,%rbx,2)
	jne	.L204
.L203:
	cmpl	$0, -92(%rsp)
	jle	.L217
	movslq	%r11d,%rcx
	leal	-1(%r11), %eax
	movq	%rdx, %r8
	salq	$2, %rcx
	xorl	%r15d, %r15d
	movl	$0, -96(%rsp)
	movq	%rcx, -80(%rsp)
	leaq	4(,%rax,4), %rbx
	leaq	-56(%rsp), %rcx
.L216:
	movq	-88(%rsp), %rdi
	movq	-72(%rsp), %r10
	movl	(%rdi,%r15), %eax
	movl	(%r10,%r15), %ebp
	leal	-7(%rax), %edx
	cmpl	%edx, %ebp
	jge	.L206
	movq	6408(%rsp), %rdx
	movslq	%ebp,%rax
	movq	-88(%rsp), %rdi
	movq	%rsi, -64(%rsp)
	leaq	(%rdx,%rax,8), %r14
	leaq	16(%rdx,%rax,8), %r13
	leaq	32(%rdx,%rax,8), %r12
	leaq	48(%rdx,%rax,8), %r10
	.p2align 4,,7
	.p2align 3
.L207:
	testl	%r11d, %r11d
	jle	.L293
	movl	(%r8), %eax
	movlpd	(%rcx), %xmm1
	leal	-1(%r11), %edx
	movl	$1, %esi
	movhpd	8(%rcx), %xmm1
	andl	$1, %edx
	addl	%ebp, %eax
	cmpl	$1, %r11d
	cltq
	movapd	%xmm1, %xmm0
	movupd	(%r9,%rax,8), %xmm2
	movupd	16(%r9,%rax,8), %xmm3
	mulpd	%xmm1, %xmm2
	movupd	32(%r9,%rax,8), %xmm4
	mulpd	%xmm1, %xmm3
	mulpd	%xmm1, %xmm4
	movupd	48(%r9,%rax,8), %xmm1
	movl	$4, %eax
	mulpd	%xmm0, %xmm1
	jle	.L210
	testl	%edx, %edx
	je	.L208
	movl	4(%r8), %eax
	movlpd	16(%rcx), %xmm9
	movl	$2, %esi
	movhpd	24(%rcx), %xmm9
	addl	%ebp, %eax
	cmpl	$2, %r11d
	cltq
	movupd	(%r9,%rax,8), %xmm8
	movupd	16(%r9,%rax,8), %xmm6
	mulpd	%xmm9, %xmm8
	movupd	32(%r9,%rax,8), %xmm7
	mulpd	%xmm9, %xmm6
	movupd	48(%r9,%rax,8), %xmm5
	movl	$8, %eax
	mulpd	%xmm9, %xmm7
	addpd	%xmm8, %xmm2
	mulpd	%xmm9, %xmm5
	addpd	%xmm6, %xmm3
	addpd	%xmm7, %xmm4
	addpd	%xmm5, %xmm1
	jle	.L210
	.p2align 4,,7
	.p2align 3
.L208:
	movl	(%r8,%rax), %edx
	movlpd	(%rcx,%rax,4), %xmm0
	addl	$2, %esi
	movlpd	16(%rcx,%rax,4), %xmm10
	movhpd	8(%rcx,%rax,4), %xmm0
	addl	%ebp, %edx
	movhpd	24(%rcx,%rax,4), %xmm10
	movslq	%edx,%rdx
	movupd	(%r9,%rdx,8), %xmm15
	movupd	16(%r9,%rdx,8), %xmm14
	mulpd	%xmm0, %xmm15
	movupd	32(%r9,%rdx,8), %xmm13
	mulpd	%xmm0, %xmm14
	movupd	48(%r9,%rdx,8), %xmm11
	movl	4(%r8,%rax), %edx
	mulpd	%xmm0, %xmm13
	addpd	%xmm15, %xmm2
	addq	$8, %rax
	mulpd	%xmm0, %xmm11
	addl	%ebp, %edx
	addpd	%xmm14, %xmm3
	movslq	%edx,%rdx
	cmpl	%esi, %r11d
	movupd	(%r9,%rdx,8), %xmm9
	addpd	%xmm13, %xmm4
	movupd	16(%r9,%rdx,8), %xmm8
	addpd	%xmm11, %xmm1
	mulpd	%xmm10, %xmm9
	movupd	32(%r9,%rdx,8), %xmm6
	mulpd	%xmm10, %xmm8
	movupd	48(%r9,%rdx,8), %xmm5
	mulpd	%xmm10, %xmm6
	addpd	%xmm9, %xmm2
	mulpd	%xmm10, %xmm5
	addpd	%xmm8, %xmm3
	addpd	%xmm6, %xmm4
	addpd	%xmm5, %xmm1
	jg	.L208
.L210:
	movupd	%xmm2, (%r14)
	addl	$8, %ebp
	addq	$64, %r14
	movupd	%xmm3, (%r13)
	addq	$64, %r13
	movupd	%xmm4, (%r12)
	addq	$64, %r12
	movupd	%xmm1, (%r10)
	addq	$64, %r10
	movl	(%rdi,%r15), %eax
	leal	-7(%rax), %esi
	cmpl	%ebp, %esi
	jg	.L207
	movq	-64(%rsp), %rsi
.L206:
	cmpl	%ebp, %eax
	jle	.L211
	movslq	%ebp,%rdx
	decl	%eax
	movq	6408(%rsp), %r13
	leaq	1(%rdx), %r10
	subl	%ebp, %eax
	xorl	%ebp, %ebp
	leaq	(%r10,%rax), %r12
	.p2align 4,,7
	.p2align 3
.L212:
	movq	%rbp, -112(%rsp)
	testl	%r11d, %r11d
	movlpd	-112(%rsp), %xmm0
	jle	.L215
	movslq	(%r8),%rax
	leaq	-4(%rbx), %rdi
	shrq	$2, %rdi
	andl	$7, %edi
	addq	%rdx, %rax
	movlpd	(%r9,%rax,8), %xmm0
	movl	$4, %eax
	cmpq	%rbx, %rax
	mulsd	(%rsi), %xmm0
	je	.L215
	testq	%rdi, %rdi
	je	.L213
	cmpq	$1, %rdi
	je	.L273
	cmpq	$2, %rdi
	.p2align 4,,5
	.p2align 3
	je	.L274
	cmpq	$3, %rdi
	.p2align 4,,5
	.p2align 3
	je	.L275
	cmpq	$4, %rdi
	.p2align 4,,5
	.p2align 3
	je	.L276
	cmpq	$5, %rdi
	.p2align 4,,5
	.p2align 3
	je	.L277
	cmpq	$6, %rdi
	.p2align 4,,5
	.p2align 3
	je	.L278
	movslq	4(%r8),%r14
	movl	$8, %eax
	addq	%rdx, %r14
	movlpd	(%r9,%r14,8), %xmm9
	mulsd	8(%rsi), %xmm9
	addsd	%xmm9, %xmm0
.L278:
	movslq	(%r8,%rax),%rdi
	addq	%rdx, %rdi
	movlpd	(%r9,%rdi,8), %xmm10
	mulsd	(%rsi,%rax,2), %xmm10
	addq	$4, %rax
	addsd	%xmm10, %xmm0
.L277:
	movslq	(%r8,%rax),%r14
	addq	%rdx, %r14
	movlpd	(%r9,%r14,8), %xmm11
	mulsd	(%rsi,%rax,2), %xmm11
	addq	$4, %rax
	addsd	%xmm11, %xmm0
.L276:
	movslq	(%r8,%rax),%rdi
	addq	%rdx, %rdi
	movlpd	(%r9,%rdi,8), %xmm12
	mulsd	(%rsi,%rax,2), %xmm12
	addq	$4, %rax
	addsd	%xmm12, %xmm0
.L275:
	movslq	(%r8,%rax),%r14
	addq	%rdx, %r14
	movlpd	(%r9,%r14,8), %xmm13
	mulsd	(%rsi,%rax,2), %xmm13
	addq	$4, %rax
	addsd	%xmm13, %xmm0
.L274:
	movslq	(%r8,%rax),%rdi
	addq	%rdx, %rdi
	movlpd	(%r9,%rdi,8), %xmm14
	mulsd	(%rsi,%rax,2), %xmm14
	addq	$4, %rax
	addsd	%xmm14, %xmm0
.L273:
	movslq	(%r8,%rax),%r14
	addq	%rdx, %r14
	movlpd	(%r9,%r14,8), %xmm15
	mulsd	(%rsi,%rax,2), %xmm15
	addq	$4, %rax
	cmpq	%rbx, %rax
	addsd	%xmm15, %xmm0
	je	.L215
.L213:
	movslq	(%r8,%rax),%rdi
	movslq	4(%r8,%rax),%r14
	addq	%rdx, %rdi
	addq	%rdx, %r14
	movlpd	(%r9,%rdi,8), %xmm8
	movslq	8(%r8,%rax),%rdi
	movlpd	(%r9,%r14,8), %xmm6
	movslq	12(%r8,%rax),%r14
	mulsd	(%rsi,%rax,2), %xmm8
	mulsd	8(%rsi,%rax,2), %xmm6
	addq	%rdx, %rdi
	movlpd	(%r9,%rdi,8), %xmm7
	addq	%rdx, %r14
	movslq	16(%r8,%rax),%rdi
	movlpd	(%r9,%r14,8), %xmm5
	movslq	20(%r8,%rax),%r14
	mulsd	16(%rsi,%rax,2), %xmm7
	addsd	%xmm8, %xmm0
	mulsd	24(%rsi,%rax,2), %xmm5
	addq	%rdx, %rdi
	movlpd	(%r9,%rdi,8), %xmm1
	addq	%rdx, %r14
	movslq	24(%r8,%rax),%rdi
	movlpd	(%r9,%r14,8), %xmm4
	movslq	28(%r8,%rax),%r14
	addsd	%xmm6, %xmm0
	mulsd	32(%rsi,%rax,2), %xmm1
	mulsd	40(%rsi,%rax,2), %xmm4
	addq	%rdx, %rdi
	movlpd	(%r9,%rdi,8), %xmm3
	addq	%rdx, %r14
	addsd	%xmm7, %xmm0
	movlpd	(%r9,%r14,8), %xmm2
	mulsd	48(%rsi,%rax,2), %xmm3
	mulsd	56(%rsi,%rax,2), %xmm2
	addq	$32, %rax
	cmpq	%rbx, %rax
	addsd	%xmm5, %xmm0
	addsd	%xmm1, %xmm0
	addsd	%xmm4, %xmm0
	addsd	%xmm3, %xmm0
	addsd	%xmm2, %xmm0
	jne	.L213
.L215:
	cmpq	%r12, %r10
	movsd	%xmm0, (%r13,%rdx,8)
	movq	%r10, %rdx
	je	.L211
	incq	%r10
	jmp	.L212
.L211:
	incl	-96(%rsp)
	addq	$4, %r15
	addq	-80(%rsp), %r8
	movl	-96(%rsp), %ebp
	cmpl	%ebp, -92(%rsp)
	jg	.L216
.L217:
	addq	$6352, %rsp
	popq	%rbx
	popq	%rbp
	popq	%r12
	popq	%r13
	popq	%r14
	popq	%r15
	ret
.L293:
	xorpd	%xmm1, %xmm1
	movapd	%xmm1, %xmm4
	movapd	%xmm1, %xmm3
	movapd	%xmm1, %xmm2
	jmp	.L210
	.cfi_endproc
.LFE527:
	.size	doperate_as_, .-doperate_as_
	.ident	"GCC: (Debian 4.4.0-5) 4.4.0"
	.section	.note.GNU-stack,"",@progbits
