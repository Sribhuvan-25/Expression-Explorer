/**
 * Publication-ready PNG export for the hand-built SVG charts. Both
 * BoxPlot and SurvivalCurve style themselves with var(--token) CSS custom
 * properties (via style="stroke:var(--rule)" and Tailwind's fill-[var(--x)]
 * classes), which only resolve inside the live DOM — a cloned/serialized
 * SVG rasterized off-screen has no access to those variables and would
 * render with black/broken paint. So this walks the *live* SVG and its
 * clone in lockstep, reading each live element's resolved fill/stroke/color
 * via getComputedStyle and writing the literal value onto the clone.
 */

interface ExportOptions {
  title: string;
  subtitle?: string;
  statLines: string[];
  scale?: number;
}

function inlineComputedPaint(liveRoot: SVGSVGElement, clonedRoot: SVGSVGElement): void {
  const walk = (live: Element, cloned: Element) => {
    if (live instanceof SVGElement) {
      const computed = getComputedStyle(live);
      // Always write the resolved value, regardless of whether it came
      // from an inline attribute, an inline style, or (as with Tailwind's
      // fill-[var(--x)] arbitrary-value classes) a stylesheet rule with no
      // matching attribute at all — the clone has none of those sources
      // once it's serialized outside the live document, so every relevant
      // property has to be captured as a literal here.
      const fill = computed.fill;
      if (fill && fill !== "none") cloned.setAttribute("fill", fill);
      const stroke = computed.stroke;
      if (stroke && stroke !== "none") cloned.setAttribute("stroke", stroke);
      if (live instanceof SVGTextElement) {
        (cloned as SVGTextElement).setAttribute("font-family", computed.fontFamily);
      }
    }
    const liveChildren = Array.from(live.children);
    const clonedChildren = Array.from(cloned.children);
    liveChildren.forEach((child, i) => {
      if (clonedChildren[i]) walk(child, clonedChildren[i]);
    });
  };
  walk(liveRoot, clonedRoot);
}

export async function exportSvgAsPng(svgEl: SVGSVGElement, filename: string, options: ExportOptions): Promise<void> {
  const scale = options.scale ?? 2;
  const viewBox = svgEl.viewBox.baseVal;
  const chartWidth = viewBox.width || svgEl.clientWidth;
  const chartHeight = viewBox.height || svgEl.clientHeight;

  const clone = svgEl.cloneNode(true) as SVGSVGElement;
  inlineComputedPaint(svgEl, clone);
  clone.setAttribute("width", String(chartWidth));
  clone.setAttribute("height", String(chartHeight));
  // Hover tooltips render as sibling HTML, not inside the SVG, but strip
  // any stray zero-opacity hover-target circles' event handlers moot since
  // this is a static clone — no cleanup needed there.

  const bodyStyle = getComputedStyle(document.body);
  const surfaceColor = bodyStyle.getPropertyValue("--surface").trim() || "#ffffff";
  const inkColor = bodyStyle.getPropertyValue("--ink").trim() || "#131A22";
  const inkMuteColor = bodyStyle.getPropertyValue("--ink-mute").trim() || "#5C6B79";
  const ruleColor = bodyStyle.getPropertyValue("--rule").trim() || "#DCE2E8";
  const fontFamily = (bodyStyle.getPropertyValue("--f-body").trim() || "sans-serif").replace(/"/g, "");
  const monoFamily = (bodyStyle.getPropertyValue("--f-mono").trim() || "monospace").replace(/"/g, "");

  const headerH = 68;
  const statsH = 20 + options.statLines.length * 18 + 12;
  const footerH = 28;
  const padding = 24;
  const totalW = chartWidth + padding * 2;
  const totalH = headerH + chartHeight + statsH + footerH + padding * 2;

  const canvas = document.createElement("canvas");
  canvas.width = totalW * scale;
  canvas.height = totalH * scale;
  const ctx = canvas.getContext("2d");
  if (!ctx) throw new Error("Canvas 2D context unavailable");
  ctx.scale(scale, scale);

  ctx.fillStyle = surfaceColor;
  ctx.fillRect(0, 0, totalW, totalH);

  ctx.fillStyle = inkColor;
  ctx.font = `700 17px ${fontFamily}`;
  ctx.fillText(options.title, padding, padding + 22);

  if (options.subtitle) {
    ctx.fillStyle = inkMuteColor;
    ctx.font = `400 11px ${fontFamily}`;
    ctx.fillText(options.subtitle, padding, padding + 40);
  }

  ctx.strokeStyle = ruleColor;
  ctx.lineWidth = 1;
  ctx.beginPath();
  ctx.moveTo(padding, headerH + padding - 8);
  ctx.lineTo(totalW - padding, headerH + padding - 8);
  ctx.stroke();

  const svgString = new XMLSerializer().serializeToString(clone);
  const svgBlob = new Blob([svgString], { type: "image/svg+xml;charset=utf-8" });
  const svgUrl = URL.createObjectURL(svgBlob);

  const img = new Image();
  img.width = chartWidth;
  img.height = chartHeight;
  try {
    await new Promise<void>((resolve, reject) => {
      img.onload = () => resolve();
      img.onerror = () => reject(new Error("Failed to rasterize chart SVG"));
      img.src = svgUrl;
    });
    ctx.drawImage(img, padding, headerH + padding, chartWidth, chartHeight);
  } finally {
    URL.revokeObjectURL(svgUrl);
  }

  let statsY = headerH + chartHeight + padding + 20;
  ctx.fillStyle = inkColor;
  ctx.font = `600 11px ${monoFamily}`;
  for (const line of options.statLines) {
    ctx.fillText(line, padding, statsY);
    statsY += 18;
  }

  ctx.strokeStyle = ruleColor;
  ctx.beginPath();
  ctx.moveTo(padding, totalH - footerH);
  ctx.lineTo(totalW - padding, totalH - footerH);
  ctx.stroke();

  ctx.fillStyle = inkMuteColor;
  ctx.font = `400 9px ${monoFamily}`;
  ctx.fillText("Expression Explorer · exported " + new Date().toISOString().slice(0, 10), padding, totalH - 10);

  const pngBlob = await new Promise<Blob | null>((resolve) => canvas.toBlob(resolve, "image/png"));
  if (!pngBlob) throw new Error("Failed to encode PNG");

  const downloadUrl = URL.createObjectURL(pngBlob);
  const a = document.createElement("a");
  a.href = downloadUrl;
  a.download = filename.endsWith(".png") ? filename : `${filename}.png`;
  document.body.appendChild(a);
  a.click();
  a.remove();
  URL.revokeObjectURL(downloadUrl);
}
