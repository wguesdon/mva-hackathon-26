-- REVTeX two-column mode rejects longtable. This filter writes each table as a plain centered
-- tabular. Quarto places a table with a cross-reference label (#tbl-...) in its own float with the
-- caption, and that float already exists when this filter runs, so the filter never adds a float.
--
-- Fitting a table into a REVTeX column needs three decisions, taken from the width of the widest cell
-- of each column, measured in characters.
--
--   1. A table that fits the column is written as it is.
--   2. A table that is too wide for a column but fits the page is wrapped in \begin{widetext}, the
--      REVTeX environment that spans both columns.
--   3. A table that is too wide for the page keeps widetext and gives each long column a \parbox, so
--      that its text wraps. A p column is not the fix: REVTeX redefines the row handling of tabular
--      and a p column fails there with "Extra \or", with or without an array package prefix.
--
-- The capacities below are characters of the body font, measured from the rendered output.
local COLUMN_CAPACITY = 44
local PAGE_CAPACITY = 92
local WRAP_AT = 20
local USABLE = 0.95

local function cell_latex(cell)
  local doc = pandoc.Pandoc(cell.contents)
  local text = pandoc.write(doc, "latex")
  return (text:gsub("\n+$", ""):gsub("\n", " "))
end

local function all_rows(tbl)
  local rows = {}
  for _, row in ipairs(tbl.head.rows) do
    table.insert(rows, row)
  end
  for _, body in ipairs(tbl.bodies) do
    for _, row in ipairs(body.body) do
      table.insert(rows, row)
    end
  end
  return rows
end

-- The longest rendered cell of each column, measured in characters.
local function column_widths(tbl)
  local widths = {}
  for index = 1, #tbl.colspecs do
    widths[index] = 0
  end
  for _, row in ipairs(all_rows(tbl)) do
    for index, cell in ipairs(row.cells) do
      local length = #cell_latex(cell)
      if widths[index] and length > widths[index] then
        widths[index] = length
      end
    end
  end
  return widths
end

-- Two characters of padding sit on each side of every column.
local function natural_width(widths)
  local total = 0
  for _, width in ipairs(widths) do
    total = total + width + 3
  end
  return total
end

-- The \parbox width of each column that must wrap, as a share of the line. The columns that fit keep
-- their natural width, and the rest of the line is shared between the long columns in proportion to
-- their content. An empty result means that no column needs a \parbox.
local function wrap_shares(widths, capacity)
  local narrow, wide = 0, 0
  for _, width in ipairs(widths) do
    if width > WRAP_AT then
      wide = wide + width
    else
      narrow = narrow + width + 3
    end
  end
  local shares = {}
  local available = USABLE - narrow / capacity
  if wide == 0 or available <= 0.1 then
    return shares
  end
  for index, width in ipairs(widths) do
    if width > WRAP_AT then
      shares[index] = available * width / wide
    end
  end
  return shares
end

local function row_latex(row, shares)
  local cells = {}
  for index, cell in ipairs(row.cells) do
    local text = cell_latex(cell)
    if shares[index] then
      text = string.format("\\parbox[t]{%.3f\\linewidth}{\\raggedright %s}", shares[index], text)
    end
    table.insert(cells, text)
  end
  return table.concat(cells, " & ") .. " \\\\"
end

local function align_spec(tbl)
  local map = { AlignLeft = "l", AlignRight = "r", AlignCenter = "c", AlignDefault = "l" }
  local spec = {}
  for _, colspec in ipairs(tbl.colspecs) do
    table.insert(spec, map[colspec[1]] or "l")
  end
  return table.concat(spec)
end

function Table(tbl)
  if not FORMAT:match("latex") then
    return nil
  end
  local widths = column_widths(tbl)
  local natural = natural_width(widths)
  local wide_page = natural > COLUMN_CAPACITY
  local shares = {}
  if natural > PAGE_CAPACITY then
    shares = wrap_shares(widths, PAGE_CAPACITY)
  end
  local lines = {}
  if wide_page then
    table.insert(lines, "\\begin{widetext}")
  end
  table.insert(lines, "\\begin{center}")
  table.insert(lines, "\\begin{tabular}{" .. align_spec(tbl) .. "}")
  table.insert(lines, "\\hline\\hline")
  for _, row in ipairs(tbl.head.rows) do
    table.insert(lines, row_latex(row, shares))
  end
  table.insert(lines, "\\hline")
  for _, body in ipairs(tbl.bodies) do
    for _, row in ipairs(body.body) do
      table.insert(lines, row_latex(row, shares))
    end
  end
  table.insert(lines, "\\hline\\hline")
  table.insert(lines, "\\end{tabular}")
  table.insert(lines, "\\end{center}")
  if wide_page then
    table.insert(lines, "\\end{widetext}")
  end
  return pandoc.RawBlock("latex", table.concat(lines, "\n"))
end
