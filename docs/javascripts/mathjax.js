// MathJax setup for the docs site: render $...$ and $$...$$ written in the Markdown files.
// (GitHub renders the same syntax by itself, so one source works in both places.)
// Reference: https://squidfunk.github.io/mkdocs-material/reference/math/
window.MathJax = {
  tex: { inlineMath: [["\\(", "\\)"]], displayMath: [["\\[", "\\]"]], processEscapes: true, processEnvironments: true },
  options: { ignoreHtmlClass: ".*|", processHtmlClass: "arithmatex" },
};
// re-typeset after Material's instant navigation swaps the page
document$.subscribe(() => {
  MathJax.startup.output.clearCache();
  MathJax.typesetClear();
  MathJax.texReset();
  MathJax.typesetPromise();
});
