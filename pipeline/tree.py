from __future__ import annotations
from .parse_file import Profile


class Node:
    children: list[Node]

    def __init__(self, profile: Profile) -> None:
        self.profile = profile

    def add_child(self, node: Node):
        self.children.append(node)


class Tree:
    root: Node
