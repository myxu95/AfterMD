#!/usr/bin/env python3
"""
测试不同的splitch链组命名模式
帮助验证AfterMD的链识别逻辑
"""

import re

def test_chain_patterns():
    """测试各种可能的链组命名模式"""
    
    print("🧪 AfterMD 链组模式识别测试")
    print("=" * 50)
    
    # 测试用例：各种可能的splitch输出格式
    test_cases = [
        # 标准格式
        ("Protein_chain_A", "字母链标识"),
        ("Protein_chain_B", "字母链标识"), 
        ("Protein_chain_C", "字母链标识"),
        
        # 数字格式
        ("Protein_chain_1", "数字链标识"),
        ("Protein_chain_2", "数字链标识"),
        ("Protein_chain_3", "数字链标识"),
        
        # 紧凑格式
        ("Protein_chainA", "紧凑字母格式"),
        ("Protein_chainB", "紧凑字母格式"),
        ("ProteinA", "极简字母格式"),
        ("ProteinB", "极简字母格式"),
        ("Protein1", "极简数字格式"),
        ("Protein2", "极简数字格式"),
        
        # ch格式
        ("ch0_Protein", "ch数字格式"),
        ("ch1_Protein", "ch数字格式"),
        ("ch2_Protein", "ch数字格式"),
        
        # Chain格式
        ("Chain_A", "Chain字母格式"),
        ("Chain_B", "Chain字母格式"), 
        ("Chain_1", "Chain数字格式"),
        ("Chain_2", "Chain数字格式"),
        
        # 应该被排除的组
        ("System", "系统组 - 应排除"),
        ("Protein", "整体蛋白质 - 应排除"),
        ("Water", "水分子 - 应排除"),
        ("SOL", "溶剂 - 应排除"),
        ("Ion", "离子 - 应排除"),
        ("Backbone", "骨架 - 应排除"),
        ("MainChain", "主链 - 应排除"),
        
        # 边界情况
        ("protein_chain_test", "包含protein和chain"),
        ("custom_protein_chain_X", "自定义含protein和chain"),
        ("some_chain_group", "只含chain"),
        ("protein_only", "只含protein"),
    ]
    
    # 定义匹配模式（与实际代码相同）
    def is_chain_group_name(group_name: str) -> bool:
        group_lower = group_name.lower()
        
        # 精确模式匹配
        chain_patterns = [
            r'^protein_chain_[a-z]$',       # Protein_chain_A
            r'^protein_chain_\d+$',         # Protein_chain_1
            r'^protein_chain[a-z]$',        # Protein_chainA
            r'^protein_chain\d+$',          # Protein_chain1
            r'^ch\d+_protein$',             # ch0_Protein
            r'^chain_[a-z]$',               # Chain_A
            r'^chain_\d+$',                 # Chain_1
            r'^protein[a-z]$',              # ProteinA
            r'^protein\d+$',                # Protein1
        ]
        
        # 检查精确模式
        for pattern in chain_patterns:
            if re.match(pattern, group_lower):
                return True, f"精确匹配: {pattern}"
        
        # 宽泛匹配：包含protein和chain
        if 'protein' in group_lower and 'chain' in group_lower:
            return True, "宽泛匹配: 包含protein和chain"
        
        # 排除明显的非链组
        non_chain_keywords = ['system', 'water', 'sol', 'ion', 'sod', 'cla', 'backbone', 'mainchain']
        for keyword in non_chain_keywords:
            if keyword in group_lower:
                return False, f"排除: 包含关键词'{keyword}'"
        
        return False, "无匹配"
    
    # 运行测试
    print("\n📋 测试结果:")
    print("-" * 80)
    print(f"{'组名':<25} {'描述':<20} {'匹配结果':<8} {'匹配原因'}")
    print("-" * 80)
    
    total_tests = len(test_cases)
    correct_matches = 0
    
    for group_name, description in test_cases:
        is_match, reason = is_chain_group_name(group_name)
        
        # 判断是否为预期结果
        expected_match = "chain" in description.lower() and "应排除" not in description
        is_correct = (is_match == expected_match)
        
        if is_correct:
            correct_matches += 1
            status = "✅"
        else:
            status = "❌"
        
        result = "是" if is_match else "否"
        print(f"{group_name:<25} {description:<20} {result:<8} {reason}")
    
    print("-" * 80)
    print(f"\n📊 测试统计:")
    print(f"总测试数: {total_tests}")
    print(f"正确匹配: {correct_matches}")
    print(f"准确率: {correct_matches/total_tests*100:.1f}%")
    
    # 显示匹配的链组
    print(f"\n🎯 识别为链组的:")
    for group_name, description in test_cases:
        is_match, reason = is_chain_group_name(group_name)
        if is_match:
            print(f"  • {group_name} ({description})")
    
    print(f"\n🚫 排除的组:")
    for group_name, description in test_cases:
        is_match, reason = is_chain_group_name(group_name)
        if not is_match:
            print(f"  • {group_name} ({description})")

def demonstrate_usage():
    """演示实际使用场景"""
    print(f"\n\n🔧 实际使用示例:")
    print("=" * 50)
    
    # 模拟不同GROMACS版本的splitch输出
    scenarios = [
        {
            "name": "GROMACS 2023 - 字母链",
            "groups": ["System", "Protein", "Protein_chain_A", "Protein_chain_B", "Water"]
        },
        {
            "name": "GROMACS 2022 - 数字链", 
            "groups": ["System", "Protein", "Protein_chain_1", "Protein_chain_2", "SOL"]
        },
        {
            "name": "GROMACS 2021 - ch格式",
            "groups": ["System", "Protein", "ch0_Protein", "ch1_Protein", "ch2_Protein", "Ion"]
        },
        {
            "name": "自定义系统",
            "groups": ["System", "ProteinA", "ProteinB", "ProteinC", "Water", "Backbone"]
        }
    ]
    
    def is_chain_simple(group_name):
        group_lower = group_name.lower()
        patterns = [r'^protein_chain_[a-z]$', r'^protein_chain_\d+$', r'^ch\d+_protein$', r'^protein[a-z]$']
        return any(re.match(p, group_lower) for p in patterns) or \
               ('protein' in group_lower and 'chain' in group_lower)
    
    for scenario in scenarios:
        print(f"\n📋 {scenario['name']}:")
        chains = [g for g in scenario['groups'] if is_chain_simple(g)]
        print(f"   发现链组: {chains}")
        if chains:
            # 模拟选择最短链（这里假设第一个最短）
            print(f"   选择最短: {chains[0]}")

if __name__ == "__main__":
    test_chain_patterns()
    demonstrate_usage()
    
    print(f"\n\n💡 使用建议:")
    print("- 如果发现新的链组命名格式，请更新匹配模式")
    print("- 在实际运行中观察日志输出，确认识别正确")
    print("- 可以通过设置DEBUG日志级别查看详细匹配过程")